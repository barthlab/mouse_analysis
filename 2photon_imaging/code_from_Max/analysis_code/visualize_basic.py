import matplotlib
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from utils import report_info
from Trials import *


color_list = ['blue', 'red', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
color_list = color_list + color_list
color_dict = {"Puff": "green", "Cued_Puff": "green", "UnCued_Puff": "green", "Blank": "yellow", "Beep": "purple",
              "No Beep": "blue",
              "ToneOn": "red", "ToneOff": "red",
              "PuffOn": "red", "PuffOff": "red",
              "NoBeepOn": "yellow", "NoBeepOff": "yellow",
              "BlankOn": "blue", "BlankOff": "blue",
              "CR+": "red", "CR-": "green", "Neurite": "gray"}
ls_dict = {"Puff": "-", "Cued_Puff": "-", "UnCued_Puff": "--", "Blank": "-", "Beep": "-", "No Beep": "-",
           "ToneOn": "-", "ToneOff": "-", "PuffOn": "-", "PuffOff": "-",
           "NoBeepOn": "-", "NoBeepOff": "-", "BlankOn": "-", "BlankOff": "-",}
lw_dict = {"CR+": 1, "CR-": 1, "Neurite": 1}
alpha_dict = {"CR+": 0.7, "CR-": 0.7, "Neurite": 0.2}

cmap = matplotlib.cm.get_cmap('plasma')
matplotlib.rcParams.update({'font.size': 11})
dy_ratio = 99.9


######################
# Basic visualization
######################


def simple_raster(raw_calcium, save_dir):
    num_cell, trace_len = raw_calcium.shape

    fig, ax = plt.subplots(2, 1, figsize=(10, 6), gridspec_kw={'height_ratios': [0.2, 1]},
                           sharex='all')
    sns.heatmap(raw_calcium, cmap='coolwarm', cbar=False, ax=ax[0], vmin=0, vmax=2)

    x = np.arange(trace_len)
    y_shift = 0
    for cell_id in range(num_cell):
        ax[1].plot(x, raw_calcium[cell_id] + y_shift, color=color_list[cell_id], lw=0.2)
        y_shift += np.percentile(raw_calcium[cell_id], dy_ratio)
    ax[1].spines[['right', 'top', 'left', 'bottom']].set_visible(False)
    ax[0].set_xticks([], [])
    ax[1].set_xticks([], [])
    fig.savefig(save_dir, bbox_inches='tight', dpi=300)
    plt.close(fig)
    report_info(f"Fig saved: {save_dir}", "b")


def extract_dy(all_data):
    array_data = np.array(all_data)
    return np.percentile(array_data-np.min(array_data), dy_ratio)[np.newaxis]


def trials_raster(trials: TrialsCollection, save_dir):
    num_session = len(set(trials.extract("session_name")))
    dy = np.concatenate([np.zeros(1),
                         np.percentile(np.array(trials.extract("df_f")), dy_ratio, axis=(0, 2)),
                         extract_dy(trials.extract("movement")),
                         # np.percentile(np.array(trials.extract("xy_off")), dy_ratio, axis=(0, 1))[np.newaxis],
                         ])
    dy = np.cumsum(dy, axis=0)
    dx = 10  # frames
    fig, ax = plt.subplots(num_session, 1, figsize=(12, 2 * num_session), sharex='all', sharey='all')
    if num_session == 1:
        ax = [ax, ]
    for session_id in range(num_session):
        tmp_trials = trials.filter(session_id=session_id)
        x_offset = 0
        session_name = "None"
        for trial in tmp_trials:
            session_name = trial.session_name
            tmp_activity, events = trial.df_f, trial.events
            sep_len = tmp_activity.shape[-1]
            x = np.arange(sep_len) + x_offset
            num_cell = tmp_activity.shape[0]
            for cell_id in range(num_cell):
                ax[session_id].plot(x, tmp_activity[cell_id] + dy[cell_id], color=color_list[cell_id], lw=0.2)
            ax[session_id].plot(x, trial.movement + dy[num_cell], color='orange', lw=0.2, alpha=0.7)
            ax[session_id].axhline(0 + dy[num_cell], xmin=0.05, xmax=0.95, color='black', lw=0.05, alpha=0.7, ls='--')
            # ax[session_id].plot(x, trial.xy_off + dy[num_cell + 1] + 1, color=color_list[-2], lw=0.2, alpha=0.7)
            for event, pos in events.items():
                ax[session_id].axvline(x=x_offset + pos, color=color_dict[event], lw=0.2, alpha=0.7)
            x_offset += sep_len + dx

        ax[session_id].spines[['right', 'top', 'left', 'bottom']].set_visible(False)
        ax[session_id].set_title(session_name)
        ax[session_id].set_xticks([], [])
    fig.savefig(save_dir, bbox_inches='tight', dpi=500)
    plt.close(fig)
    report_info(f"Fig saved: {save_dir}", "b")


def session_raster(sessions: TrialsCollection, save_dir):

    lw_mod = 2
    beh_off = 1.2
    num_session = len(sessions)
    # dy = np.concatenate([np.zeros(1),
    #                      np.percentile(np.array(sessions.extract("df_f")), dy_ratio, axis=(0, 2))*1.2+1.,
    #                      extract_dy(sessions.extract("movement")),
    #                      # np.percentile(np.array(sessions.extract("xy_off")), dy_ratio, axis=(0, 1))[np.newaxis],
    #                      ])
    # dy = np.clip(dy, a_min=0, a_max=10)
    # dy[-2] += 2
    # dy = np.cumsum(dy, axis=0)
    fig, ax = plt.subplots(num_session, 1, figsize=(8, 2 * num_session), sharex='all', sharey='all')
    if num_session == 1:
        ax = [ax, ]
    for session_id in range(num_session):
        tmp_session = sessions[session_id]

        ####
        dy = np.concatenate([np.zeros(1),
                             np.percentile(np.array(tmp_session.get("df_f")), dy_ratio, axis=1) * 1.2 + 1.,
                             extract_dy(tmp_session.get("movement")),
                             # np.percentile(np.array(sessions.extract("xy_off")), dy_ratio, axis=(0, 1))[np.newaxis],
                             ])
        dy = np.clip(dy, a_min=0, a_max=10)
        dy[-2] += 2
        dy = np.cumsum(dy, axis=0)

        ####


        tmp_activity, trials_type, trials_frame = tmp_session.df_f, tmp_session.trials_type, tmp_session.trials_frame
        x = np.arange(tmp_activity.shape[-1])
        num_cell = tmp_activity.shape[0]
        for cell_id in range(num_cell):
            ax[session_id].plot(x, tmp_activity[cell_id] + dy[cell_id], color=color_list[cell_id], lw=0.1*lw_mod)
        ax[session_id].plot(x, tmp_session.movement + dy[num_cell], color='blue', lw=0.1*lw_mod, alpha=0.7)
        ax[session_id].axhline(dy[num_cell], xmin=0.05, xmax=0.95, color='black', lw=0.05*lw_mod, alpha=0.7, ls='--')

        ax[session_id].scatter(x, tmp_session.position + dy[num_cell] + beh_off, facecolor='purple', edgecolor='none', s=0.1)
        # ax[session_id].axhline(dy[num_cell] + beh_off, xmin=0.05, xmax=0.95, color='black', lw=0.05*lw_mod, alpha=0.7, ls='--')
        # ax[session_id].axhline(dy[num_cell] + beh_off+1, xmin=0.05, xmax=0.95, color='black', lw=0.05*lw_mod, alpha=0.7, ls='--')

        # ax[session_id].plot(x, tmp_session.grooming + dy[num_cell] + 2*beh_off, color='green', lw=0.1*lw_mod)
        #
        # ax[session_id].scatter(x, tmp_session.whisking + dy[num_cell] + 3*beh_off, facecolor='red', edgecolor='none', s=0.1)
        # ax[session_id].axhline(dy[num_cell] + 3*beh_off, xmin=0.05, xmax=0.95, color='black', lw=0.05*lw_mod, alpha=0.7, ls='--')
        # ax[session_id].axhline(dy[num_cell] + 3*beh_off+1, xmin=0.05, xmax=0.95, color='black', lw=0.05*lw_mod, alpha=0.7, ls='--')

        for trial_id, x_offset in enumerate(trials_frame):
            # ax[session_id].axvline(x=x_offset, ymin=0., ymax=1, color=color_dict[trials_type[trial_id]], lw=0.01*lw_mod)
            # ax[session_id].axvline(x=x_offset, ymin=0.95, ymax=1, color=color_dict[trials_type[trial_id]], lw=0.2*lw_mod)
            ax[session_id].axvline(x=x_offset, ymin=0., ymax=1, color='gray', lw=0.01*lw_mod)
            ax[session_id].axvline(x=x_offset, ymin=0.98, ymax=1, color='gray', lw=0.2*lw_mod)
        ax[session_id].spines[['right', 'top', 'left', 'bottom']].set_visible(False)
        ax[session_id].set_title(tmp_session.session_name)
        ax[session_id].set_xticks([], [])
    fig.savefig(save_dir, bbox_inches='tight', dpi=500)
    plt.close(fig)
    report_info(f"Fig saved: {save_dir}", "b")
