import matplotlib
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from utils import report_info, read_video_dir, read_xlsx, read_csv_dir
from Trials import *
import os
import os.path as path
import cv2

tmp_dir = path.join(".", "_tmp")
os.makedirs(tmp_dir, exist_ok=True)

px = 1 / plt.rcParams['figure.dpi']

debug = False
video_fps = 30
session_time = 6e5
color_list = ['blue', 'red', 'green', 'yellow', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
color_list = color_list + color_list
color_dict = {"Puff": "green", "Cued_Puff": "green", "UnCued_Puff": "green", "Blank": "yellow", "Beep": "purple",

              "No Beep": "blue",
              "ToneOn": "red", "ToneOff": "red",
              "PuffOn": "red", "PuffOff": "red",
              "NoBeepOn": "yellow", "NoBeepOff": "yellow",
              "BlankOn": "blue", "BlankOff": "blue",
              "CR+": "red", "CR-": "green", "Neurite": "gray"}


def activity2img(session: Trial, start_id: int, end_id: int, w, h, minor_marker=True):
    random_name_token = np.random.randint(100)
    tmp_activity, trials_type, trials_frame = session.df_f, session.trials_type, session.trials_frame
    dy = np.concatenate([np.zeros(1), np.percentile(tmp_activity, 99.9, axis=1) * 1.1, ])
    dy = np.cumsum(np.clip(dy, a_min=0, a_max=10), axis=0)
    fig, ax = plt.subplots(1, 1, figsize=(w * px, h * px), facecolor='black')
    end_id = tmp_activity.shape[-1] if end_id == -1 else end_id
    x = np.arange(start_id, end_id)
    num_cell = tmp_activity.shape[0]
    for cell_id in range(num_cell):
        ax.plot(x, tmp_activity[cell_id][start_id: end_id] + dy[cell_id],
                color='white', lw=1)
        ax.axhline(dy[cell_id], xmin=0., xmax=1., color='gray', lw=0.5, alpha=0.7)

    for trial_id, x_offset in enumerate(trials_frame):
        if start_id <= x_offset < end_id:
            if minor_marker:
                ax.axvline(x=x_offset, ymin=0., ymax=1, color=color_dict[trials_type[trial_id]], lw=0.7)
            else:
                ax.axvline(x=x_offset, ymin=0.95, ymax=1, color=color_dict[trials_type[trial_id]], lw=0.7)
    # ax.spines[['right', 'top', 'left', 'bottom']].set_visible(False)
    ax.set_axis_off()
    ax.set_xlim(start_id, end_id)
    ax.set_ylim(-1, dy[-1])
    fig.savefig(path.join(tmp_dir, f"{random_name_token}.jpg"), bbox_inches='tight',
                pad_inches=0)  # pad_inches = 0 important
    plt.close(fig)
    image = cv2.imread(path.join(tmp_dir, f"{random_name_token}.jpg"))
    return image


def behavior2img(session: Trial, start_id: int, end_id: int, w, h, minor_marker=True):
    random_name_token = np.random.randint(100)
    tmp_activity, trials_type, trials_frame = session.df_f, session.trials_type, session.trials_frame
    fig, ax = plt.subplots(1, 1, figsize=(w * px, h * px), facecolor='black')
    end_id = len(session.movement) if end_id == -1 else end_id
    x = np.arange(start_id, end_id)

    ax.scatter(x, session.movement[start_id: end_id],
               facecolor='yellow', edgecolor='none', s=5)
    ax.axhline(0, xmin=0., xmax=1., color='white', lw=0.5, alpha=0.7)

    ax.scatter(x, session.position[start_id: end_id] + 2,
               facecolor='cyan', edgecolor='none', s=5)
    ax.axhline(2, xmin=0., xmax=1., color='white', lw=0.5, alpha=0.7)
    # ax.axhline(3, xmin=0., xmax=1., color='white', lw=0.5, alpha=0.7)

    ax.scatter(x, session.grooming[start_id: end_id] + 3.5,
               facecolor='green', edgecolor='none', s=5)
    ax.axhline(3.5, xmin=0., xmax=1., color='white', lw=0.5, alpha=0.7)
    # ax.axhline(4.5, xmin=0., xmax=1., color='white', lw=0.5, alpha=0.7)

    ax.scatter(x, session.whisking[start_id: end_id] + 5,
               facecolor='red', edgecolor='none', s=5)
    ax.axhline(5, xmin=0., xmax=1., color='white', lw=0.5, alpha=0.7)
    # ax.axhline(6, xmin=0., xmax=1., color='white', lw=0.5, alpha=0.7)

    for trial_id, x_offset in enumerate(trials_frame):
        if start_id <= x_offset < end_id:
            if minor_marker:
                ax.axvline(x=x_offset, ymin=0., ymax=1, color=color_dict[trials_type[trial_id]], lw=0.7)
            else:
                ax.axvline(x=x_offset, ymin=0.95, ymax=1, color=color_dict[trials_type[trial_id]], lw=0.7)
    # ax.spines[['right', 'top', 'left', 'bottom']].set_visible(False)
    ax.set_axis_off()
    ax.set_xlim(start_id, end_id)
    ax.set_ylim(np.min(session.movement) - 0.5, 6.5)
    fig.savefig(path.join(tmp_dir, f"{random_name_token}.jpg"), bbox_inches='tight',
                pad_inches=0)  # pad_inches = 0 important
    plt.close(fig)
    image = cv2.imread(path.join(tmp_dir, f"{random_name_token}.jpg"))
    return image


def highlight_image(ori_image, start_id, end_id, frame_data):
    new_image = (ori_image + 0.).astype(float)
    height, width, _ = new_image.shape
    new_image[:, :int(width * start_id / frame_data)] *= 0.5
    new_image[:, int(width * end_id / frame_data):] *= 0.5
    return new_image.astype(np.uint8)


def render_frames(tmp_session, tmp_video, tmp_offset, half_window):
    num_frame_data = tmp_session.df_f.shape[-1]
    num_frame_video, video_height, video_width = tmp_video.shape
    global_activity = activity2img(tmp_session, start_id=0, end_id=-1, w=1000, h=250, minor_marker=False)
    global_behavior = behavior2img(tmp_session, start_id=0, end_id=-1, w=1000, h=250, minor_marker=False)
    sidebar_height, sidebar_width, _ = global_behavior.shape

    video2data = (np.linspace(tmp_offset, tmp_offset + 1000 * num_frame_video / video_fps,
                              num_frame_video) * num_frame_data / session_time).astype(int)

    final_video_height = int(video_height * sidebar_width / video_width)
    annotation_side = np.zeros((final_video_height, sidebar_width, 3), np.uint8)
    unit_height = int(sidebar_height / 5)
    unit_width = int(sidebar_width/3)
    cv2.putText(annotation_side, "whisking", (unit_width, unit_height), cv2.FONT_HERSHEY_SIMPLEX,
                1, (0, 0, 255), 2)
    cv2.putText(annotation_side, "grooming", (unit_width, 2 * unit_height), cv2.FONT_HERSHEY_SIMPLEX,
                1, (0, 255, 0), 2)
    cv2.putText(annotation_side, "position", (unit_width, 3 * unit_height), cv2.FONT_HERSHEY_SIMPLEX,
                1, (255, 255, 0), 2)
    cv2.putText(annotation_side, "velocity", (unit_width, 4 * unit_height), cv2.FONT_HERSHEY_SIMPLEX,
                1, (0, 255, 255), 2)
    cv2.putText(annotation_side, "Events", (unit_width, 7 * unit_height), cv2.FONT_HERSHEY_SIMPLEX,
                1, (255, 255, 255), 2)
    cv2.putText(annotation_side, "Vertical Puff", (unit_width, 8 * unit_height), cv2.FONT_HERSHEY_SIMPLEX,
                1, (0, 0, 255), 2)
    cv2.putText(annotation_side, "Horizontal Puff", (unit_width, 9 * unit_height), cv2.FONT_HERSHEY_SIMPLEX,
                1, (255, 0, 0), 2)
    cv2.putText(annotation_side, "Cell activity", (unit_width, 12 * unit_height), cv2.FONT_HERSHEY_SIMPLEX,
                1, (255, 255, 255), 2)

    for cur_video_frame_id, cur_data_frame_id in enumerate(video2data):
        if cur_data_frame_id <= half_window or cur_data_frame_id + half_window >= num_frame_data:
            continue
        window_dict = {"start_id": cur_data_frame_id - half_window, "end_id": cur_data_frame_id + half_window}
        cur_video_frame = cv2.flip(cv2.resize(cv2.cvtColor(tmp_video[cur_video_frame_id], cv2.COLOR_GRAY2BGR),
                                              None, fx=sidebar_width / video_width, fy=sidebar_width / video_width), 0)
        local_activity = activity2img(tmp_session, **window_dict, w=1000, h=250, minor_marker=True)
        local_behavior = behavior2img(tmp_session, **window_dict, w=1000, h=250, minor_marker=True)

        cur_frame = np.concatenate([
            local_behavior,
            cur_video_frame,
            local_activity,
        ], axis=0)
        cur_height, cur_width, _ = cur_frame.shape
        cv2.drawMarker(cur_frame, (int(cur_width / 2), sidebar_height),
                       color=(125, 125, 125), thickness=2, markerType=cv2.MARKER_TRIANGLE_UP)
        cv2.drawMarker(cur_frame, (int(cur_width / 2), cur_height - sidebar_height),
                       color=(125, 125, 125), thickness=2, markerType=cv2.MARKER_TRIANGLE_DOWN)

        cur_annotation = np.concatenate([
            highlight_image(global_behavior, **window_dict, frame_data=num_frame_data),
            annotation_side,
            highlight_image(global_activity, **window_dict, frame_data=num_frame_data),
        ], axis=0)
        final_frame = np.concatenate([cur_annotation, cur_frame], axis=1)
        yield final_frame


def tmp_session_videos_interface(sessions: TrialsCollection, mice_dir, save_dir, half_window: int = 100):
    os.makedirs(save_dir, exist_ok=True)
    video_dir = path.join(mice_dir, 'videos')
    rpi_time = read_xlsx(path.join(mice_dir, f"time_point.xlsx"), 0)
    ttl = read_csv_dir(path.join(mice_dir, f"ttl"), 0)
    all_videos = read_video_dir(video_dir)
    num_session = len(sessions)

    time_offset = [rpi_time[session_id][0, 0] - (ttl[session_id][1, 1] - ttl[session_id][0, 1])
                   for session_id in range(num_session)]
    for session_id in range(num_session):
        fourcc = cv2.VideoWriter_fourcc(*'mp4v')
        video_save_path = path.join(save_dir, all_videos[session_id].split('mouse_video_')[-1].replace("npy", "mp4"))
        out = cv2.VideoWriter(video_save_path, fourcc, 30.0, (1550, 935))
        frame_cnt = 0
        for frame in render_frames(sessions[session_id], np.load(all_videos[session_id]), time_offset[session_id],
                                   half_window):
            print(f"Frame cnt: {frame_cnt}")
            if debug:
                print(frame.shape)
                cv2.imshow('frame', frame)
                cv2.waitKey(1)
            out.write(frame)
            frame_cnt += 1
        out.release()
        cv2.destroyAllWindows()
        if debug:
            exit()
