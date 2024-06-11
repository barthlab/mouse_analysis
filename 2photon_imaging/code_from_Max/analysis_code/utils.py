from scipy.io import loadmat, savemat
import numpy as np
import pandas as pd
import os
import os.path as path
from scipy.stats import ttest_ind
from sklearn.feature_selection import f_regression
from sklearn.linear_model import LinearRegression, Lasso, HuberRegressor
from colorist import Color
import cv2


def z_score_np(a, axis):
    return (a - np.mean(a, axis=axis, keepdims=True)) / np.std(a, axis=axis, keepdims=True)


# def sem(a, axis=None):
#     tmp_a = np.array(a)
#     if axis is None:
#         return np.std(tmp_a)
#     else:
#         return np.std(tmp_a, axis=axis)


def sem(a, axis=None):
    tmp_a = np.array(a)
    if axis is None:
        return np.std(tmp_a)/np.sqrt(len(tmp_a.flatten()))
    else:
        return np.std(tmp_a, axis=axis)/np.sqrt(tmp_a.shape[axis])


def elaborate_mat(mat_array: dict):
    keys_list = list(mat_array.keys())
    print(f"Containing:")
    for key in keys_list:
        if isinstance(mat_array[key], np.ndarray):
            additional_info = str(mat_array[key].shape)
        elif isinstance(mat_array[key], str):
            additional_info = mat_array[key]
        elif isinstance(mat_array[key], list):
            additional_info = f"{len(mat_array[key])} elements: F10 {str(mat_array[key][:10])}"
        elif isinstance(mat_array[key], bytes):
            additional_info = str(mat_array[key])
        else:
            additional_info = ''
        print(f"{key}: {type(mat_array[key])} : {additional_info}")


def load_mat(mat_dir, show_details=False):
    mat_array = loadmat(mat_dir)
    print(f"\n{mat_dir} Loaded")
    if show_details:
        elaborate_mat(mat_array)
    return mat_array


def save_mat(mat_array, save_dir: str, show_details=False):
    savemat(save_dir, mat_array)
    print(f"\n{save_dir} Saved")
    if show_details:
        elaborate_mat(mat_array)


def read_xlsx(table_dir, header=0, show_details=False, output_sheet_name=False):
    xl_file = pd.ExcelFile(table_dir)
    xlsx_dict = {sheet_id: xl_file.parse(sheet_name, header=header).to_numpy()
                 for sheet_id, sheet_name in enumerate(xl_file.sheet_names)}
    if show_details:
        elaborate_mat(xlsx_dict)
    if output_sheet_name:
        sheet_names = {sheet_id: sheet_name for sheet_id, sheet_name in enumerate(xl_file.sheet_names)}
        return xlsx_dict, sheet_names
    else:
        return xlsx_dict


def read_csv_dir(table_dir, header=0, show_details=False, output_sheet_name=False):
    raw_dict = {}
    for (dirpath, dirnames, filenames) in os.walk(table_dir):
        for filename in filenames:
            csv_file = pd.read_csv(path.join(table_dir, filename))
            raw_dict[filename] = csv_file.to_numpy()
    sheet_names = {sheet_id: sheet_name for sheet_id, sheet_name in enumerate(sorted(raw_dict.keys()))}
    xlsx_dict = {sheet_id: raw_dict[sheet_names[sheet_id]] for sheet_id in sheet_names.keys()}
    if show_details:
        elaborate_mat(xlsx_dict)
    if output_sheet_name:
        return xlsx_dict, sheet_names
    else:
        return xlsx_dict


def read_npy_dir(npy_dir, show_details=False, output_sheet_name=False):
    raw_dict = {}
    for (dirpath, dirnames, filenames) in os.walk(npy_dir):
        for filename in filenames:
            raw_dict[filename] = np.load(path.join(npy_dir, filename))
    sheet_names = {sheet_id: sheet_name for sheet_id, sheet_name in enumerate(sorted(raw_dict.keys()))}
    npy_dict = {sheet_id: raw_dict[sheet_names[sheet_id]] for sheet_id in sheet_names.keys()}
    if show_details:
        elaborate_mat(npy_dict)
    if output_sheet_name:
        return npy_dict, sheet_names
    else:
        return npy_dict


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


def read_video_dir(video_dir, show_details=False):
    video_list = []
    for root, dirs, files in os.walk(video_dir):
        for name in files:
            if name.split(".")[-1] == "h264":
                video_list.append(os.path.join(root, name))
    if show_details:
        print(video_list)
    all_videos = []
    for video_path in video_list:
        if not os.path.exists(video_path.replace("h264", "npy")):
            print(f"Retrieving Frames from {video_path}")
            cap = cv2.VideoCapture(video_path)
            all_frames = []
            frame_cnt = 0
            while True:
                ret0, frame = process_capture(cap)
                if not ret0:
                    break
                if show_details:
                    print(f"{frame_cnt} frames")
                all_frames.append(frame)
                frame_cnt += 1
            all_frames = np.array(all_frames)
            np.save(video_path.replace("h264", "npy"), all_frames)
        all_videos.append(video_path.replace("h264", "npy"))
    return all_videos


def p_value(v1, v2):
    assert v1.ndim == 1 == v2.ndim
    return ttest_ind(v1, v2).pvalue


def ftest(xs, ys):
    assert len(xs) == len(ys)
    statis, p = f_regression(xs[:, np.newaxis], ys)
    reg = HuberRegressor().fit(xs[:, np.newaxis], ys)
    # reg = LinearRegression().fit(xs[:, np.newaxis], ys)
    return reg.intercept_, float(reg.coef_), float(p)


def pearson_cor(x: np.ndarray, y: np.ndarray):
    assert x.shape == y.shape
    x_hat = x - np.mean(x)
    y_hat = y - np.mean(y)
    x_hat = x_hat / np.sqrt(np.sum(x_hat ** 2))
    y_hat = y_hat / np.sqrt(np.sum(y_hat ** 2))
    return np.sum(x_hat * y_hat)


def huber_cor(x: np.ndarray, y: np.ndarray):
    assert x.shape == y.shape
    reg = HuberRegressor().fit(x.flatten()[:, np.newaxis] / np.std(x), y.flatten() / np.std(y))
    return float(reg.coef_)


def get_id_ood_indices(x_vector, y_vector, thrd):
    ood_indices = np.argwhere(np.logical_or(np.abs(x_vector - np.mean(x_vector)) > thrd * np.std(x_vector),
                                            np.abs(y_vector - np.mean(y_vector)) > thrd * np.std(y_vector)))[:, 0]
    id_indices = np.argwhere(np.logical_and(np.abs(x_vector - np.mean(x_vector)) <= thrd * np.std(x_vector),
                                            np.abs(y_vector - np.mean(y_vector)) <= thrd * np.std(y_vector)))[:, 0]
    return ood_indices, id_indices


def arg_sorted_where(names, key: str, day_labels=None):
    if day_labels is None:
        day_labels = ["ACC1)", "ACC2)", "ACC3)", "ACC4)", "ACC5)", "ACC6)", ]
        # "SAT1)", "SAT2)", "SAT3)", "SAT4)", "SAT5)", "SAT6)", "SAT7)", "SAT8)", "SAT9)", "SAT10)"]
    selected_names = []
    for day in day_labels:
        for name in names:
            if (key in name) and (day in name):
                selected_names.append(name)
                break
    indices = [names.index(name) for name in selected_names]
    return selected_names, indices



def report_info(info, color="red"):
    if color == "red" or color == "r":
        print(f"{Color.RED}{info}{Color.OFF}")
    elif color == "blue" or color == "b":
        print(f"{Color.BLUE}{info}{Color.OFF}")
    elif color == "yellow" or color == "y":
        print(f"{Color.YELLOW}{info}{Color.OFF}")
    elif color == "green" or color == "g":
        print(f"{Color.GREEN}{info}{Color.OFF}")
    elif color == "cyan" or color == "c":
        print(f"{Color.CYAN}{info}{Color.OFF}")
    else:
        print(f"{Color.MAGENTA}{info}{Color.OFF}")


def simple_beeswarm2(y, nbins=None, width=1.):
    """
    Returns x coordinates for the points in ``y``, so that plotting ``x`` and
    ``y`` results in a bee swarm plot.
    """
    y = np.asarray(y)
    if nbins is None:
        # nbins = len(y) // 6
        nbins = np.ceil(len(y) / 0.01).astype(int)

    # Get upper bounds of bins
    x = np.zeros(len(y))

    nn, ybins = np.histogram(y, bins=nbins)
    nmax = nn.max()

    #Divide indices into bins
    ibs = []#np.nonzero((y>=ybins[0])*(y<=ybins[1]))[0]]
    for ymin, ymax in zip(ybins[:-1], ybins[1:]):
        i = np.nonzero((y>ymin)*(y<=ymax))[0]
        ibs.append(i)

    # Assign x indices
    dx = width / (nmax // 2)
    for i in ibs:
        yy = y[i]
        if len(i) > 1:
            j = len(i) % 2
            i = i[np.argsort(yy)]
            a = i[j::2]
            b = i[j+1::2]
            x[a] = (0.5 + j / 3 + np.arange(len(b))) * dx
            x[b] = (0.5 + j / 3 + np.arange(len(b))) * -dx

    return x