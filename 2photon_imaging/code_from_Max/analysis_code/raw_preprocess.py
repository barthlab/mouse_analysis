import numpy as np

from utils import *
from visualize_basic import simple_raster, trials_raster, session_raster
from Trials import *
from scipy.ndimage import gaussian_filter

trial_range = (25, 25)
baseline_range = 15
session_time = 6e5
unit_frame_speed = 5.11 * 29.8 / 6000  # dm/s
tpr = 600


def drop_frames(mat_array, offset_threshold=20, block_len_threshold=2):
    ops = mat_array['ops'][0][0]
    F = mat_array['F']
    Fneu = mat_array['Fneu']
    num_neurites, total_frames = F.shape
    xy_off = np.sqrt(ops['xoff'] ** 2 + ops['yoff'] ** 2)[0]
    xy_off = (np.abs(xy_off - np.mean(xy_off)) / 10).astype(np.float32)

    x_off = np.argwhere(np.abs(ops['xoff']) > offset_threshold)[:, 1]
    y_off = np.argwhere(np.abs(ops['yoff']) > offset_threshold)[:, 1]
    corrupted_frames = np.unique(np.concatenate([x_off, y_off]))
    corrupted_frames = np.concatenate([corrupted_frames, [3 * total_frames]])
    block_start = 0
    dropped_frames = np.zeros(total_frames)
    for i in range(len(corrupted_frames) - 1):
        cur_frame, nxt_frame = corrupted_frames[i], corrupted_frames[i + 1]
        if nxt_frame > cur_frame + 1:
            block_len = i - block_start + 1
            start_frame = corrupted_frames[block_start]

            if (block_len <= block_len_threshold) and (start_frame + block_len < total_frames):
                interp = np.linspace(0, 1, block_len + 2)[np.newaxis, 1:-1]
                F[:, start_frame:start_frame + block_len] = (interp * F[:, start_frame + block_len, np.newaxis] +
                                                             (1 - interp) * F[:, start_frame - 1, np.newaxis])
                Fneu[:, start_frame:start_frame + block_len] = (interp * Fneu[:, start_frame + block_len, np.newaxis] +
                                                                (1 - interp) * Fneu[:, start_frame - 1, np.newaxis])
            else:
                dropped_frames[start_frame:start_frame + block_len] = 1
            block_start = i + 1

    mat_array['F'], mat_array['Fneu'] = F, Fneu
    return mat_array, dropped_frames, xy_off


def df_f_calculation(mat_array):
    F = mat_array['F']
    Fneu = mat_array['Fneu']
    is_cell = np.argwhere(mat_array['iscell'][:, 0] == 1)[:, 0]
    print(f"Cell Index : {is_cell}")
    central_activity, background_activity = F[is_cell], Fneu[is_cell]
    calibrated_activity = central_activity - 0.7 * background_activity
    return calibrated_activity


def trial_extraction(calcium, time_point_time, ttl_time, dropped_frames, xy_off, session_name, num_session,
                     behavior: dict[np.ndarray],
                     additional_param_names: tuple = (), global_baseline=True, **kwargs) -> (TrialsCollection, TrialsCollection):
    """

    :param global_baseline:
    :param additional_param_names:
    :param calcium: shape (num_cell, total_frames)
    :param time_point_time: extracted from time_point file, first column corresponding to time points received in microscope
    :param ttl_time: extracted from ttl file, first column represent to trial type, second column represent trial time generated in raspberry Pi
    :param dropped_frames: shape (total_frames) 1 for drop, 0 for reserve
    :param xy_off: shape (total_frames)
    :param session_name: shape (num_session, #params) manual parameter for each session
    :param behavior: dictionary of different behavior data, (num_session, #session frames) each
    :param num_session:
    :return:
    """
    # prepare
    all_trial, all_session = TrialsCollection(), TrialsCollection()
    num_cell, num_total_frames = calcium.shape
    num_session_frames = int(num_total_frames / num_session)
    num_add_params = len(additional_param_names)

    # reshape
    session_calcium = np.reshape(calcium, (num_cell, num_session_frames, num_session), order='F')
    session_drop_frames = np.reshape(dropped_frames, (num_session_frames, num_session), order='F')
    session_xy_off = np.reshape(xy_off, (num_session_frames, num_session), order='F')

    # main extraction loop
    for session_id in range(num_session):
        # align ttl time for each trial to time point time
        trials_time = ttl_time[session_id][1:, 1] - ttl_time[session_id][1, 1] + time_point_time[session_id][0, 0]
        precise_frame = (num_session_frames * trials_time / session_time) - 1  # minus 1 for signal delay
        # convert to frame
        trials_frame = np.floor(precise_frame).astype(int)

        # extract other parameters
        trials_type = ttl_time[session_id][1:, 0]
        trials_add_params = ttl_time[session_id][1:, 2:]
        tmp_calcium = session_calcium[..., session_id]
        tmp_behavior = {key: behavior[key][session_id] for key in behavior.keys()}
        num_trials = -1
        tmp_baseline = []

        for trial_index, trial_f in enumerate(trials_frame):
            # extract trial info
            flu = tmp_calcium[:, trial_f - trial_range[0]: trial_f + trial_range[1]]
            baseline = np.mean(flu[:, trial_range[0] - baseline_range:trial_range[0]], axis=-1, keepdims=True)
            behavior_trial = {key: tmp_behavior[key][trial_f - trial_range[0]: trial_f + trial_range[1]]
                              for key in tmp_behavior.keys()}

            tmp_xy_off = session_xy_off[trial_f - trial_range[0]: trial_f + trial_range[1], session_id]

            tmp_baseline.append(baseline)
            if "On" in trials_type[trial_index]:
                num_trials += 1
                drop_flag = (np.sum(session_drop_frames[trial_f - trial_range[0]:trial_f + trial_range[1],
                                    session_id]) > 0)
                tmp_events = {trials_type[trial_index]: trial_range[0] + precise_frame[trial_index] - trial_f,
                              trials_type[trial_index + 1]: trial_range[0] + precise_frame[trial_index + 1] - trial_f}
                if not global_baseline:
                    flu = (flu-baseline)/baseline

                new_trial = Trial(df_f=flu, trial_type=trials_type[trial_index][:-2], drop=drop_flag,
                                  events=tmp_events, trial_id=num_trials, xy_off=tmp_xy_off,
                                  exp_mice=session_name[session_id, 0], exp_day=session_name[session_id, 1],
                                  fov=session_name[session_id, 2], info=session_name[session_id, 3],
                                  session_id=session_id,
                                  session_name=f" | ".join(map(str, list(session_name[session_id]) + [session_id, ])),
                                  **behavior_trial,
                                  **{param_name: trials_add_params[trial_index, param_id]
                                     for param_id, param_name in enumerate(additional_param_names)})
                all_trial.add(new_trial)

        # normalize to baseline
        tmp_baseline = np.mean(np.stack(tmp_baseline, axis=0), axis=0)
        extract_calcium = (tmp_calcium - tmp_baseline) / tmp_baseline
        if global_baseline:
            for trial in all_trial.filter(exp_mice=session_name[session_id, 0], exp_day=session_name[session_id, 1],
                                          fov=session_name[session_id, 2], session_id=session_id).src_data:
                trial.update(df_f=(trial.df_f - tmp_baseline) / tmp_baseline)
        new_session = Trial(df_f=extract_calcium, trials_type=trials_type, trials_frame=precise_frame,
                            xy_off=session_xy_off[..., session_id],
                            exp_mice=session_name[session_id, 0], exp_day=session_name[session_id, 1],
                            fov=session_name[session_id, 2], info=session_name[session_id, 3],
                            session_id=session_id,
                            session_name=f" | ".join(map(str, list(session_name[session_id]) + [session_id, ])),
                            **tmp_behavior,
                            **{param_name: trials_add_params[:, param_id]
                               for param_id, param_name in enumerate(additional_param_names)})
        all_session.add(new_session)
    return all_trial, all_session


def movement_extraction(distance_dict, num_session, session_frames, time_offset):
    movement_dict, position_dict = {}, {}
    for session_id in range(num_session):
        time_points = distance_dict[session_id][:, 0]
        positions = distance_dict[session_id][:, 1]
        directions = distance_dict[session_id][:, 2]
        tmp_movement = np.zeros(session_frames, dtype=np.float32)
        tmp_position = np.zeros(session_frames, dtype=np.float32)
        tmp_points = np.zeros(session_frames, dtype=np.float32)
        time_frames = np.floor(session_frames * (time_points - time_offset[session_id]) / session_time).astype(int)

        for record_id in range(1, len(time_points)):
            if session_frames > time_frames[record_id] >= 0:
                tmp_points[time_frames[record_id]] += 1
                tmp_position[time_frames[record_id]] += positions[record_id]
        tmp_position = tmp_position/np.clip(tmp_points, a_min=1, a_max=None)

        cur_pos = tmp_position[np.argwhere(tmp_points > 0)[0, 0]]
        for frame_id in range(session_frames):
            if tmp_points[frame_id] == 0 and tmp_position[frame_id] == 0:
                tmp_position[frame_id] = cur_pos
            elif tmp_points[frame_id] > 0:
                cur_pos = tmp_position[frame_id]
            else:
                raise NotImplementedError
        tmp_position = gaussian_filter(tmp_position, sigma=1)
        tmp_movement[1:] = tmp_position[1:] - tmp_position[:-1]
        movement_dict[session_id] = (tmp_movement * unit_frame_speed).astype(np.float32)
        position_dict[session_id] = (tmp_position % tpr)/tpr
    return movement_dict, position_dict


def behavior_extraction(behavior_dict, num_session, session_frames, time_offset):
    video_fps = 30
    new_behavior_dict = {}
    for session_id in range(num_session):
        raw_trace = behavior_dict[session_id]
        xp = np.linspace(time_offset[session_id], time_offset[session_id]+1000*len(raw_trace)/video_fps, len(raw_trace))
        x = np.linspace(0, session_time, session_frames)
        new_behavior_dict[session_id] = np.interp(x, xp, raw_trace)
    return new_behavior_dict


def data_fetcher(mice_dir, num_session, save_dir, elaborate: bool = True, **kwargs) -> (TrialsCollection, TrialsCollection):
    print(f"Fetching {mice_dir}")
    os.makedirs(save_dir, exist_ok=True)
    raw_data = load_mat(path.join(mice_dir, 'Fall.mat'), show_details=elaborate)
    rpi_time = read_xlsx(path.join(mice_dir, f"time_point.xlsx"), 0, elaborate)
    ttl = read_csv_dir(path.join(mice_dir, f"ttl"), 0, elaborate)
    session_name = read_xlsx(path.join(mice_dir, f"session_names.xlsx"), None, elaborate)[0]
    raw_mov = read_csv_dir(path.join(mice_dir, "distance"), 0, elaborate)
    # raw_whisking = read_npy_dir(path.join(mice_dir, "whisking"), elaborate)
    # raw_grooming = read_npy_dir(path.join(mice_dir, "grooming"), elaborate)

    raw_data, dropped_pos, xyoff = drop_frames(raw_data)
    calcium = df_f_calculation(raw_data)

    ttl_rpi_offset = [ttl[session_id][1, 1] - rpi_time[session_id][0, 0] for session_id in range(num_session)]
    video_rpi_offset = [rpi_time[session_id][0, 0] - (ttl[session_id][1, 1]-ttl[session_id][0, 1])
                        for session_id in range(num_session)]

    session_frames = int(calcium.shape[1] / num_session)
    movement, position = movement_extraction(raw_mov, num_session, session_frames, ttl_rpi_offset)
    # whisking = behavior_extraction(raw_whisking, num_session, session_frames, video_rpi_offset)
    # grooming = behavior_extraction(raw_grooming, num_session, session_frames, video_rpi_offset)

    trials, sessions = trial_extraction(calcium, rpi_time, ttl, dropped_pos, xyoff, session_name, num_session,
                                        behavior={ "movement": movement,"position": position,},
                                        # behavior={"movement": movement, "position": position,
                                        #           "whisking": whisking, "grooming": grooming, },
                                        **kwargs)
    if elaborate:
        simple_raster(np.concatenate(sessions.extract("df_f"), axis=-1), path.join(save_dir, "full_raster.jpg"))
        session_raster(sessions, path.join(save_dir, "session_raster.pdf"))
        trials_raster(trials, path.join(save_dir, "trials_raster.jpg"))
    return trials, sessions


def extract_peak_before(trial: Trial):
    assert trial_range[0] > 3 * baseline_range
    return np.max(trial.df_f[:, :trial_range[0] - 3 * baseline_range], axis=-1)


def extract_peak_sensitive_period(trial: Trial):
    assert trial_range[0] > baseline_range
    return np.max(trial.df_f[:, trial_range[0] - baseline_range + 1: trial_range[0]], axis=-1)


def extract_peak_first_peak(trial: Trial):
    return np.max(trial.df_f[:, trial_range[0]: trial_range[0] + 3 * baseline_range], axis=-1)


def extract_peak_second_peak(trial: Trial):
    assert trial_range[0] > 3 * baseline_range
    return np.max(trial.df_f[:, trial_range[0] + baseline_range: trial_range[0] + 4 * baseline_range], axis=-1)


def extract_peak_after(trial: Trial):
    assert trial_range[0] > 3 * baseline_range
    return np.max(trial.df_f[:, trial_range[0] + 4 * baseline_range:], axis=-1)


peak_extractor = {
    "peak_before": extract_peak_before,
    "peak_sensitive": extract_peak_sensitive_period,
    "peak_first": extract_peak_first_peak,
    "peak_second": extract_peak_second_peak,
    "peak_after": extract_peak_after
}


def extract_movement_auc(trial: Trial):
    return np.sum(np.abs(trial.movement[trial_range[0]-5: trial_range[0]+5]))


def extract_whisking_auc(trial: Trial):
    return np.sum(np.abs(trial.whisking[trial_range[0]-5: trial_range[0]+5]))


def extract_grooming_auc(trial: Trial):
    return np.sum(np.abs(trial.grooming[trial_range[0]-5: trial_range[0]+5]))


def extract_df_f_auc(trial: Trial):
    return np.sum(trial.df_f[:, trial_range[0]-5: trial_range[0]+5], axis=1)


def extract_behavior_type(trial: Trial):
    if trial.grooming_auc >=2:
        return "grooming"
    if trial.whisking_auc <=4:
        return "stationary"
    elif trial.movement_auc <= 1:
        return "whisking only"
    elif trial.whisking_auc <= 9.5:
        return "whisking while running"
    else:
        return "intense whisking"


auc_extractor = {
    "whisking_auc": extract_whisking_auc,
    "movement_auc": extract_movement_auc,
    "grooming_auc": extract_grooming_auc,
    "df_f_auc": extract_df_f_auc,
}


behavior_extractor = {
    "behavior_type": extract_behavior_type,
}