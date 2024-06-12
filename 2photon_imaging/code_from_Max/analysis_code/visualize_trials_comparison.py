import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from utils import sem, report_info
from Trials import *
import os
import os.path as path

color_list = ['blue', 'red', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
color_list = color_list + color_list
color_dict = {"Puff": "green", "Cued_Puff": "green", "UnCued_Puff": "green", "Blank": "yellow", "Beep": "purple",
              "No Beep": "blue",
              "ToneOn": "red", "ToneOff": "red", "PuffOn": "green", "PuffOff": "green",
              "NoBeepOn": "yellow", "NoBeepOff": "yellow", "BlankOn": "blue", "BlankOff": "blue",
              "CR+": "red", "CR-": "green", "Neurite": "gray"}
ls_dict = {"Puff": "-", "Cued_Puff": "-", "UnCued_Puff": "--", "Blank": "-", "Beep": "-", "No Beep": "-",
           "ToneOn": "-", "ToneOff": "-", "PuffOn": "-", "PuffOff": "-",
           "NoBeepOn": "-", "NoBeepOff": "-", "BlankOn": "-", "BlankOff": "-", }
lw_dict = {"CR+": 1, "CR-": 1, "Neurite": 1}
alpha_dict = {"CR+": 0.7, "CR-": 0.7, "Neurite": 0.2}
dy_ratio = 99.9


def mean_and_sem_plot(ax, x, trace, offset, color, lw):
    mean_trace, sem_trace = np.mean(trace, axis=0), sem(trace, axis=0)
    ax.plot(x, mean_trace + offset, color=color, lw=lw)
    ax.fill_between(x, mean_trace - sem_trace + offset, mean_trace + sem_trace + offset, color=color, lw=lw / 4,
                    alpha=0.2)


def mean_and_sem_scatter(ax, x, trace, offset, color, s):
    mean_trace, sem_trace = np.mean(trace, axis=0), sem(trace, axis=0)
    # ax.plot(x, mean_trace + offset, color=color, lw=lw)
    ax.scatter(x, mean_trace + offset, facecolor=color, edgecolor='none', s=s)

    ax.axhline(offset, xmin=0.05, xmax=0.95, color='black', lw=s/2, alpha=0.7, ls='--')
    # ax.axhline(offset + 1, xmin=0.05, xmax=0.95, color='black', lw=0.05 * lw_mod, alpha=0.7,
    #                        ls='--')
    # ax.fill_between(x, mean_trace - sem_trace + offset, mean_trace + sem_trace + offset, color=color, lw=lw / 4,
    #                 alpha=0.2)


def trials_plot(trials_axes: list[list[TrialsCollection]] | list[TrialsCollection] | TrialsCollection,
                row_names, col_names, save_dir):
    lw_mod = 4
    if isinstance(trials_axes, TrialsCollection):
        trials_axes = [[trials_axes, ], ]
    elif isinstance(trials_axes[0], TrialsCollection):
        trials_axes = [trials_axes, ]
    num_row = len(trials_axes)
    num_col = np.max([len(trials_row) for trials_row in trials_axes])
    all_together = sum(sum(trials_axes, []))
    dy = np.concatenate([np.zeros(1),
                         np.percentile(np.array(all_together.extract("df_f")), dy_ratio, axis=(0, 2)),
                         np.percentile(np.array(all_together.extract("movement"))-
                                       np.array(all_together.extract("movement")).min(),
                                       dy_ratio, axis=(0, 1))[np.newaxis],
                         np.percentile(np.array(all_together.extract("whisking")), dy_ratio, axis=(0, 1))[np.newaxis],
                         # np.percentile(np.array(all_together.extract("xy_off")), dy_ratio, axis=(0, 1))[np.newaxis],
                         ])
    dy = np.cumsum(dy, axis=0)
    tmp_alpha = .1

    fig, axs = plt.subplots(num_row, num_col, figsize=(3 * num_col, 4 * num_row), sharex='all', sharey='all')
    for row_id in range(num_row):
        for col_id in range(num_col):
            ax = axs[row_id] if num_row > 1 else axs
            ax = ax[col_id] if num_col > 1 else ax

            # ax.set_ylabel(row_names[row_id])
            ax.spines[['right', 'top', 'left', 'bottom']].set_visible(False)
            if col_id >= len(trials_axes[row_id]):
                continue
            # ax.set_xlabel(col_names[col_id])
            tmp_trials = trials_axes[row_id][col_id]
            if len(tmp_trials) == 0:
                continue

            total_df_f = np.array(tmp_trials.extract("df_f"))
            total_mov = np.array(tmp_trials.extract("movement"))
            total_whi = np.array(tmp_trials.extract("whisking"))
            # total_xy_off = np.array(tmp_trials.extract("xy_off"))
            _, num_cell, seq_len = total_df_f.shape
            x = np.arange(seq_len)

            for cell_id in range(num_cell):
                mean_and_sem_plot(ax, x, total_df_f[:, cell_id], offset=dy[cell_id], color=color_list[cell_id], lw=0.8)
            mean_and_sem_scatter(ax, x, total_mov, offset=dy[num_cell], color='blue', s=0.5)
            mean_and_sem_scatter(ax, x, total_whi, offset=dy[num_cell] + 1.2, color='red', s=0.5)
            # ax.axhline(0 + dy[num_cell], xmin=0.05, xmax=0.95, color='black', lw=0.2, alpha=0.7, ls='--')

            for trial in tmp_trials:
                tmp_df_f = trial.df_f
                for cell_id in range(num_cell):
                    ax.plot(x, tmp_df_f[cell_id] + dy[cell_id], color=color_list[cell_id],
                            lw=0.1 * lw_mod, alpha=tmp_alpha)
                # ax.plot(x, trial.movement + dy[num_cell], color='orange', lw=0.1 * lw_mod, alpha=tmp_alpha)
                # ax.plot(x, trial.xy_off + dy[num_cell + 1], color=color_list[-2], lw=0.1,
                #                         alpha=tmp_alpha)

            for event, pos in tmp_trials[0].events.items():
                ax.axvline(x=pos, color=color_dict[event], ls=ls_dict[event], lw=0.2 * lw_mod, alpha=0.7)
            # ax.set_xticks([0, 25, 50], ["-5", "0", "5"])
            ax.text(1, 1, f"n={len(tmp_trials)}", ha='right', va='top', transform=ax.transAxes)

    fig.savefig(save_dir, bbox_inches='tight', dpi=500)
    plt.close(fig)
    report_info(f"Fig saved: {save_dir}", "b")
