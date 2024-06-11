#!/bin/env python3

import argparse
import csv
import sys
import time
import os
import constants
import matplotlib.pyplot as plt
import numpy as np
import cv2
import tkinter as tk
from tkinter import filedialog
from sklearn.linear_model import HuberRegressor

# Constants for circle detection
debug = True
clicked_point = None
infer_window = 1000


def on_mouse_click(event, x, y, flags, param):
    global clicked_point
    if event == cv2.EVENT_LBUTTONDOWN:
        clicked_point = (x, y)


def process_capture(capture):
    # Grab the frame
    ret0, frame = capture.read()

    # If frame successfully grabbed
    if ret0:
        # Process frame into grayscale
        frame = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)

    return ret0, frame


def uu(array: np.ndarray):
    return array.astype(np.uint8)


def ii(array: np.ndarray):
    return array.astype(int)


def save_to_csv(radius_values_path, radius_list, center_list):
    with open(radius_values_path, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["time", "radius", "center x", "center y"])

        for idx in range(len(radius_list)):
            writer.writerow([*radius_list[idx], *center_list[idx]])


def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    """

    return np.isnan(y), lambda z: z.nonzero()[0]


def retrieve_eyes(video_path):
    cap = cv2.VideoCapture(video_path)
    cv2.namedWindow("Video")
    # cv2.setMouseCallback("Video", on_mouse_click)
    center_positions = []
    frame_cnt = 0
    while True:
        ret0, frame = process_capture(cap)
        if not ret0 or (infer_window != -1 and frame_cnt > infer_window):
            break
        print(f"{frame_cnt} frames")
        frame = uu((ii(frame) - np.min(ii(frame))) / (np.max(ii(frame)) - np.min(ii(frame))) * 255)
        screen_frame = cv2.cvtColor(frame, cv2.COLOR_GRAY2RGB)

        f_height, f_width = frame.shape[:2]
        assert (f_height, f_width) == (768, 1024)

        _, thresh_light = cv2.threshold(frame, 240, 255, cv2.THRESH_BINARY)
        contours, _ = cv2.findContours(thresh_light, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        satisfied_contours = []
        for contour in contours:
            area = cv2.contourArea(contour)
            rect = cv2.boundingRect(contour)
            x, y, width, height = rect
            radius = 0.25 * (width + height)

            symmetry_condition = (abs(1 - float(width) / float(height)) <= 0.2)
            fill_condition = (abs(1 - (area / (np.pi * radius ** 2))) <= 0.2)
            area_condition = (width * height < 12 * 12)
            # print(x, y, width, height, area, symmetry_condition, fill_condition, area_condition)

            cx, cy = int(x + width / 2), int(y + height / 2)
            pupil_crop = 10
            extract_rect_s = uu(ii(frame)[cy - pupil_crop:cy + pupil_crop, cx - pupil_crop:cx + pupil_crop])
            darkness_condition = (np.sum(ii(extract_rect_s) < 120) / (4 * pupil_crop * pupil_crop)) > 0.1
            small_condition = (width < 5 and height < 5)
            not_too_small_condition = (width >= 2 and height >= 2)
            if area_condition and darkness_condition and not_too_small_condition and(
                    fill_condition + symmetry_condition + small_condition >= 1):
                satisfied_contours.append(contour)
                center_positions.append((frame_cnt, cx, cy))
                screen_frame = cv2.rectangle(screen_frame, (cx - pupil_crop, cy - pupil_crop),
                                             (cx + pupil_crop, cy + pupil_crop), (0, 0, 255), 2)


        tmp_screen = screen_frame.copy()
        cv2.drawContours(tmp_screen, satisfied_contours, -1, (255, 0, 0), 2)
        cv2.imshow("Video", tmp_screen)
        cv2.waitKey(1)

        frame_cnt += 1
    x_pos, y_pos = HuberRegressor(epsilon=1.01), HuberRegressor(epsilon=1.01)
    center_positions = np.array(center_positions)
    x_pos.fit(center_positions[:, 0:1], center_positions[:, 1])
    y_pos.fit(center_positions[:, 0:1], center_positions[:, 2])
    print(x_pos.coef_, x_pos.intercept_)
    print(y_pos.coef_, y_pos.intercept_)
    cxs = x_pos.predict(np.arange(20000).reshape(-1, 1))
    cys = y_pos.predict(np.arange(20000).reshape(-1, 1))
    return cxs, cys


def extract_pupil(video_path, cxs, cys):
    frame_cnt = 0
    pupil_crop = 10
    eye_crop_x, eye_crop_y = 50, 35
    search_dots_num = 5

    eye_cap = cv2.VideoCapture(video_path)
    cv2.namedWindow("Eye")

    X = np.arange(eye_crop_x*2)
    Y = np.arange(eye_crop_y*2)
    X, Y = np.meshgrid(X, Y)
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator
    radius_list = []
    center_list = []
    while True:
        ret0, frame = process_capture(eye_cap)
        if not ret0 or (infer_window != -1 and frame_cnt > infer_window):
            break

        print(f"{frame_cnt} frames")
        cx, cy = int(cxs[frame_cnt]), int(cys[frame_cnt])
        pupil = uu(ii(frame)[cy - pupil_crop:cy + pupil_crop, cx - pupil_crop:cx + pupil_crop])
        eye = uu(ii(frame)[cy - eye_crop_y:cy + eye_crop_y, cx - eye_crop_x:cx + eye_crop_x])

        darkest_pts = np.array(np.unravel_index(np.argsort(ii(pupil), axis=None), pupil.shape))
        darkest_pts[0] += eye_crop_y - pupil_crop
        darkest_pts[1] += eye_crop_x - pupil_crop

        max_radius = 0
        final_radius_result = []
        final_center_result = []
        for dot_id in range(search_dots_num):
            chosen_pt = np.array([darkest_pts[0][dot_id], darkest_pts[1][dot_id]])

            # _, light_dot = cv2.threshold(pupil, 240, 255, cv2.THRESH_BINARY)
            # light_thres = np.zeros((eye_crop_y*2, eye_crop_x*2), np.uint8)
            # light_thres[eye_crop_y-pupil_crop: eye_crop_y+pupil_crop, eye_crop_x-pupil_crop: eye_crop_x+pupil_crop] = light_dot

            flood_results = []
            circle_result = []
            radius_result = []
            center_result = []
            for updiff in range(20):
                mask = np.zeros((eye_crop_y*2 + 2, eye_crop_x*2 + 2), np.uint8)
                flood_eye = eye.copy()
                cv2.floodFill(flood_eye, mask, (chosen_pt[1], chosen_pt[0]), 255, loDiff=4, upDiff=updiff)
                mask_thres = uu(np.clip(mask[1:-1, 1:-1]*255, a_max=255, a_min=0))
                # mask_thres = uu(np.clip(mask[1:-1, 1:-1]*255 + ii(light_thres), a_max=255, a_min=0))

                contours, _ = cv2.findContours(mask_thres, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
                contours = np.concatenate(contours)
                (enclosure_x, enclosure_y), enclosure_radius = cv2.minEnclosingCircle(contours)
                if enclosure_radius < 25:
                    circle_eye = eye.copy()
                    cv2.circle(circle_eye, (int(enclosure_x), int(enclosure_y)), int(enclosure_radius), (255, 0, 0), 2)
                    flood_results.append(mask_thres)
                    circle_result.append(circle_eye)
                    radius_result.append(int(enclosure_radius))
                    center_result.append((int(enclosure_x), int(enclosure_y)))
                else:
                    break
            col1 = np.concatenate([eye, *flood_results])
            col2 = np.concatenate([eye, *circle_result])
            cv2.imshow("Eye", np.concatenate([col1, col2], axis=1))
            if debug:
                cv2.waitKey(0)
            else:
                cv2.waitKey(1)

            if len(radius_result) > 0 and radius_result[-1] > max_radius:
                max_radius = radius_result[-1]
                final_radius_result = radius_result
                final_center_result = center_result

        if len(final_radius_result) > 0:
            radius_list.append(final_radius_result[-1])
            center_list.append(final_center_result[-1])
            print(final_radius_result)

            if debug:
                fig, ax = plt.subplots(1, 1, figsize=(10, 10), subplot_kw={"projection": "3d"})
                surf = ax.plot_surface(Y, X, ii(eye), cmap=cm.coolwarm, linewidth=0, antialiased=False, alpha=0.6)
                ax.scatter(darkest_pts[0][:search_dots_num], darkest_pts[1][:search_dots_num],
                           ii(eye)[darkest_pts[0][:search_dots_num], darkest_pts[1][:search_dots_num]], s=30, color='black')
                ax.zaxis.set_major_locator(LinearLocator(10))
                fig.colorbar(surf, shrink=0.5, aspect=5)
                plt.show()
        else:
            radius_list.append(np.nan)
            center_list.append((np.nan, np.nan))

        frame_cnt += 1
    return radius_list, center_list


def main(video_path, radius_path):
    cxs, cys = retrieve_eyes(video_path)
    radius_list, center_list = extract_pupil(video_path, cxs, cys)
    frame_cnt = len(radius_list)
    print(f"Total Frame Count: {frame_cnt}")
    for radius_id in range(1, frame_cnt-1):
        if (abs(radius_list[radius_id] - radius_list[radius_id-1]) > 15) and (abs(radius_list[radius_id] - radius_list[radius_id+1]) > 15):
            radius_list[radius_id] = np.nan
            center_list[radius_id] = (np.nan, np.nan)
    radius_list = np.array(radius_list)
    # nans, x = nan_helper(radius_list)
    # radius_list[nans] = np.interp(x(nans), x(~nans), radius_list[~nans])

    if frame_cnt/600 > 28:
        fps = 30
    elif frame_cnt/600 < 20:
        fps = 15
    else:
        fps = 25
    time_radius = [(frame_id/fps, radius) for frame_id, radius in enumerate(radius_list)]
    save_to_csv(radius_path, time_radius, center_list)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Detect a mouse's pupil radius given video of a mouse")
    parser.add_argument("--video_path", type=str, default="./videos", help="Path to video files")
    args = parser.parse_args()

    for root, dirs, files in os.walk(args.video_path):
        for name in files:
            if name.split(".")[-1] == "h264":
                video_dir = os.path.join(root, name)
                radius_dir = video_dir.replace(constants.VIDEO_PREFIX, constants.PUPIL_PREFIX).replace(".h264", ".csv")
                print(video_dir)
                main(video_dir, radius_dir)

