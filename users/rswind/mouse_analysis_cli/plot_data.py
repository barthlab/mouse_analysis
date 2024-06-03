import numpy as np
import seaborn as sns #v0.13.0
import matplotlib as mpl #v3.8.1
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.dates as mdates

from matplotlib.patches import Rectangle
from matplotlib.collections import LineCollection


def style_axes_helper(ax, title, xlabel, ylabel, xlim, ylim, xmajlocator, xminlocator,
               ymajlocator, yminlocator):
    ax.set_title(f"{title}", fontsize=20)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)    
    if xmajlocator != None: ax.xaxis.set_major_locator(xmajlocator)
    if xminlocator != None: ax.xaxis.set_minor_locator(xminlocator)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=18)
    if ymajlocator != None: ax.yaxis.set_major_locator(ymajlocator)
    if yminlocator != None: ax.yaxis.set_minor_locator(yminlocator)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=18)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

def style_axes(g, xlabel, ylabel, xlim, ylim, xmajmult, xminmult, ymajmult,
               yminmult, hspace=0.45, wspace=0.25, suptitle=None,title=None):
    xmaj = ticker.MultipleLocator(xmajmult) if xmajmult != None else xmajmult
    xmin = ticker.MultipleLocator(xminmult) if xminmult != None else xminmult
    ymaj = ticker.MultipleLocator(ymajmult) if ymajmult != None else ymajmult
    ymin = ticker.MultipleLocator(yminmult) if yminmult != None else yminmult

    t = title
    if g.axes_dict: 
        for (name,ax) in g.axes_dict.items():
            if title == None:
                t=name
            style_axes_helper(ax, t, xlabel, ylabel, xlim, ylim, xmaj, xmin, ymaj, ymin)
    else:
        style_axes_helper(g.figure.axes[0], '', xlabel, ylabel, xlim, ylim, xmaj, xmin, ymaj, ymin)

    g.figure.suptitle(suptitle, fontsize=25, y=1.02, x=0.5)
    g.figure.subplots_adjust(hspace=hspace, wspace=wspace)

def style_hr(g, xlabel, ylabel, xlim, ylim, hspace=0.45, wspace=0.25, 
             suptitle=None, title=None, ymajmult=2, yminmult=1, xmajmult=12):
    style_axes(g, xlabel, ylabel, xlim, ylim, xmajmult, 4, ymajmult, yminmult, 
               hspace=hspace, wspace=wspace, suptitle=suptitle, title=title)
    
    # if g.axes_dict:
    #     for (name,ax) in g.axes_dict.items():
    #         ax.add_patch(Rectangle((0, ylim[0]), xlim[1], ylim[1] - ylim[0], color="#CEE1F3", alpha=0.8, 
    #                             zorder=0, fill=True))        
    # else:
    #     ax = g.figure.axes[0]
    #     ax.add_patch(Rectangle((0, ylim[0]), xlim[1], ylim[1] - ylim[0], color="#CEE1F3", alpha=0.8, 
    #                             zorder=0, fill=True))     

def style_ms(g, xlabel, ylabel, xlim, ylim, xmajmult=500, xminmult=100, 
             ymajmult=2, yminmult=1,hspace=0.45, wspace=0.25, 
             suptitle=None,style=True):
    style_axes(g, xlabel, ylabel, xlim, ylim, xmajmult, xminmult, ymajmult, yminmult, hspace=hspace,
               wspace=wspace, suptitle=suptitle)
    if style:
        if g.axes_dict:
            for (name,ax) in g.axes_dict.items():
                ax.add_patch(Rectangle((0, ylim[0]), 500, ylim[1] - ylim[0], color="lightgrey", alpha=0.4, 
                                    zorder=0, fill=True))        
                ax.add_patch(Rectangle((1005, ylim[0]), 100, ylim[1] - ylim[0], color="#BBE9FD", alpha=0.8, 
                                        zorder=0, fill=True))
        else:
            ax = g.figure.axes[0]
            ax.add_patch(Rectangle((0, ylim[0]), 500, ylim[1] - ylim[0], color="lightgrey", alpha=0.4, 
                            zorder=0, fill=True))        
            ax.add_patch(Rectangle((1005, ylim[0]), 100, ylim[1] - ylim[0], color="#BBE9FD", alpha=0.8, 
                                    zorder=0, fill=True))


     
def style_performance(g):
    if g.axes_dict:
        for (name,ax) in g.axes_dict.items(): 
            ax.axhline(y=0, xmin=0, xmax=1, ls="-", color="black", zorder=10, lw=1)
    else:
        g.figure.axes[0].axhline(y=0, xmin=0, xmax=1, ls="-", color="black", zorder=10, lw=1)

def style_trial(g, color='k', alpha=0.3):
    if g.axes_dict: 
        for (name,ax) in g.axes_dict.items():
            ax.fill_between(ax.lines[0].get_data()[0], ax.lines[0].get_data()[1], 
                        color=color, alpha=alpha)
    else:
        ax = g.figure.axes[0]
        ax.fill_between(ax.lines[0].get_data()[0], ax.lines[0].get_data()[1], 
                        color=color, alpha=alpha)
  
def lineplot(data, x=None, y=None, errorbar=None, aspect=1.5, marker='o', 
             mec=None, ms=10, **kwargs):
    
    g = sns.relplot(data, x=x, y=y, kind="line",  aspect=aspect, 
                    errorbar=errorbar, marker=marker, mec=mec,ms=ms,
                    facet_kws={"sharey":False, "sharex":False}, **kwargs)
    return g

def barplot(data, x=None, y=None, hue=None, col=None, color=None, aspect=1,
            palette=None, hue_order=None, legend=None, col_wrap=None, 
            errorbar="se", dodge=False, capsize=0.2, sharex=False, sharey=False,
            native_scale=True, err_kws={"zorder":0.5, "lw":1}, **kwargs):
     
    g = sns.catplot(data, x=x, y=y, palette=palette, col_wrap=col_wrap,
                    hue_order=hue_order, hue=hue, color=color, col=col, 
                    kind="bar", errorbar=errorbar, dodge=dodge, capsize=capsize, 
                    err_kws=err_kws, sharex=sharex, sharey=sharey,aspect=aspect,
                    native_scale=native_scale, legend=legend, **kwargs)
    return g

def plot_hr(data, x="Time (hr)", y="lick", color='k',
            xlim=[-48, 24], ylim=[-10, 10], ymajmult=2, yminmult=1,
            ylabel="", xlabel="Time (hr)", hspace=0.45, 
            wspace=0.25, suptitle='', title=None, **kwargs):
    
    g = lineplot(data, x=x, y=y,color=color, **kwargs)
    
    style_hr(g, xlabel, ylabel, xlim, ylim, hspace=hspace, wspace=wspace, 
             suptitle=suptitle, title=title,ymajmult=ymajmult, yminmult=yminmult,)  
    
    return g 

def plot_ant_lickfreq(data, x="Time (hr)", y="lick", hue="stimulus", 
                      col="condition", palette=["green", "red"], 
                      hue_order=["stimulus", "blank"],
                      ylim=[0, 12], ylabel="Lick Frequency (Hz)", 
                      **kwargs):
    g = plot_hr(data, x=x, y=y, hue=hue, col=col, palette=palette, 
                hue_order=hue_order, ylim=ylim, ylabel=ylabel, **kwargs)
    
    return g
 
def plot_ant_perf(data, x="Time (hr)", y="lick", col="condition", **kwargs):
    g = plot_hr(data, x=x, y=y, col=col, **kwargs) 
    style_performance(g)
    return g

def plot_ant_perf_bar(data, x="Time (hr)", y="lick", hue=None, col="condition", 
                      xlim=[-50, 26], ylim=[-10, 10], palette=None, 
                      hue_order=None, color='k', ylabel="Performance", 
                      xlabel="Time (hr)", hspace=0.45, wspace=0.35, legend=True,
                      col_wrap=None, suptitle='', **kwargs):
    
    g = barplot(data, x=x, y=y, col=col, col_wrap=col_wrap, hue=hue, 
                palette=palette, hue_order=hue_order, color=color, 
                legend=legend, **kwargs)  
    style_hr(g, xlabel, ylabel, xlim, ylim, hspace=hspace, wspace=wspace, 
             suptitle=suptitle)
    style_performance(g)

    return g

def plot_lickfreq(data, x="Time (ms)", y="lick", hue="stimulus", 
                  col="Time (hr)",color=None,legend=True, col_wrap=6,
                  palette=["green", "red"], hue_order=["stimulus", "blank"], 
                  xlim=[-200, 2000], ylim=[0, 12], ylabel="Lick Frequency",
                  xlabel="Time (ms)", hspace=0.45, wspace=0.25, suptitle=None,
                  xmajmult=500, ymajmult=2, xminmult=100, 
                  yminmult=1, style=True,**kwargs):
    
    g = lineplot(data, x=x, y=y, hue=hue, col_wrap=col_wrap, col=col,
                 palette=palette, color=color, hue_order=hue_order, 
                 legend=legend, **kwargs)
    style_ms(g, xlabel, ylabel, xlim, ylim, xmajmult=xmajmult, ymajmult=ymajmult,
             xminmult=xminmult, yminmult=yminmult,hspace=hspace, wspace=wspace, 
             suptitle=suptitle,style=style)
    
    return g

def plot_perf(data, x="Time (ms)", y="lick", hue="Time (hr)", 
              col="condition", color='k', legend=True, col_wrap=None,
              palette=None, hue_order=None, xlim=[-200, 2000], ylim=[-10, 10],
              ylabel="Performance", xlabel="Time (ms)", hspace=0.45, 
              wspace=0.25, suptitle=None, **kwargs):
    g = lineplot(data, x=x, y=y, hue=hue, col_wrap=col_wrap, col=col,
                 palette=palette, color=color, hue_order=hue_order, 
                 legend=legend, **kwargs)
    style_ms(g, xlabel, ylabel, xlim, ylim, hspace=hspace, wspace=wspace, 
             suptitle=suptitle)
    style_performance(g)
    
    return g

def plot_bar(data, x=None, y=None, aspect=1, col=None, 
                   color="#3B3838", palette=None, col_wrap=None, hue_order=None,
                   hue=None, xlim=[-0.5, 2.5], ylim=[0, 600], xlabel=None, hspace=0.45, wspace=0.25,
                   ylabel="", width=0.4, ymajmult=100, yminmult=50,title='', suptitle='',
                   **kwargs):
    g = barplot(data, x=x, y=y, aspect=aspect, col=col, color=color, 
                palette=palette, hue_order=hue_order, hue=hue, 
                col_wrap=col_wrap, width=width, **kwargs)
    style_axes(g, xlabel, ylabel, xlim, ylim, None, None, ymajmult=ymajmult, 
               yminmult=yminmult, hspace=hspace, wspace=wspace, title=title, suptitle=suptitle)

    return g

def plot_trial_hr(data, x="Time (hr)", y="trial no", col="condition", color='k', 
                  hue=None, palette=None, col_wrap=None, hue_order=None, 
                  errorbar='se', err_style='bars', xlabel="Time (hr)", 
                  ylabel="Number of trials", xlim=[-48, 48], ylim=[0, 300], 
                  hspace=0.45, wspace=0.25, suptitle=None, fill=True, title=None, 
                  ymajmult=100, yminmult=50,**kwargs):
    g = lineplot(data, x=x, y=y, hue=hue, col=col, color=color, palette=palette, 
                 hue_order=hue_order, errorbar=errorbar, err_style=err_style, 
                 **kwargs)
    style_hr(g, xlabel, ylabel, xlim, ylim, ymajmult=ymajmult, yminmult=yminmult,
             hspace=hspace, wspace=wspace, suptitle=suptitle, title=title)
    # for name, ax in g.axes_dict.items():
    #     ax.spines[['right']].set_visible(True)
    #     ax.spines[['left']].set_visible(False)
    #     ax.xaxis.set_visible(False)
    #     ax.yaxis.tick_right()

    if fill: style_trial(g)
    return g

def connect_lines(ax,data, x=None,y=None, hue=None, order=None, hue_order=None):
    segs = []
    if hue and hue != x:
        for i in range(len(order)):
            hfilt = data[data[x] ==  order[i]]
            for a in hfilt["animal"].unique():
                anfilt = hfilt[hfilt["animal"] == a]
                pts = []
                for j in range(len(hue_order)):
                    op = list(anfilt[(anfilt[hue] == hue_order[j])][y])
                    if (len(op) != 0):
                        pts.append([i + 1/(len(hue_order)+1)*(j-(len(hue_order)-1)/2), op[0]])
                if(len(pts) != 0):
                
                    segs.append(pts)
    else:
        for a in data["animal"].unique():
            anfilt = data[data["animal"] == a]
            pts = []
            for j in range(len(order)):
                op = list(anfilt[(anfilt[x] == order[j])][y])
                if (len(op) != 0):
                        pts.append([j, op[0]])
            if(len(pts) != 0):
                segs.append(pts)
    linesegs = LineCollection(segs,color='k', lw=0.5, zorder=1)
    ax.add_collection(linesegs)


def plot_bar_strip(data, x=None, y=None, row=None, col=None, col_order=None, connect=False,
                   errorbar='se', color=None, pt_color=None, err_kws={'lw':1, 'zorder':0}, 
                   capsize=0.2, row_order=None,ylim=[0,10], xlim=[-0.5, 3.5], order=None, dodge=False, jitter=False,
                   hue=None, palette=None, hue_order=None, ylabel='', xlabel='',xlab_kws={}, **kwargs):
    g = sns.catplot(data, x=x, y=y, row=row,col=col,col_order=col_order,order=order,
                    errorbar=errorbar, color=color, err_kws=err_kws, capsize=capsize,kind='bar', 
                    row_order=row_order, hue=hue, palette=palette, hue_order=hue_order, **kwargs)
    if g.axes_dict:
        for name, ax in g.axes_dict.items():
            if type(name) == tuple:
                i, j = name
                cdfilt = data[(data[col] == j) & (data[row]==i)] 
                ax.set_title(f'{i} {j}')
            else:
                if col:
                    cdfilt = data[data[col] == name]
                elif row:
                    cdfilt = data[data[row] == name]
                ax.set_title(name)
            sns.stripplot(cdfilt, ax=ax, x=x, color=pt_color, y=y, 
                          dodge=dodge, jitter=jitter, hue=hue, legend=False,
                          hue_order=hue_order, palette=palette)
            if connect:
                connect_lines(ax,cdfilt, hue=hue, order=order, x=x, y=y,hue_order=hue_order)
                            
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.tick_params(labelbottom=True, labeltop=False)
            ax.set_xticklabels(ax.get_xticklabels(), **xlab_kws)
    else:
        ax = g.axes[0][0]
        sns.stripplot(data, ax=ax, x=x, color=pt_color, y=y, 
                        dodge=dodge, jitter=jitter, hue=hue, legend=False,
                        hue_order=hue_order, palette=palette)
        if connect:
            connect_lines(ax,data, hue=hue, order=order, x=x,y=y, hue_order=hue_order)
                        
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.tick_params(labelbottom=True, labeltop=False)
        ax.set_xticklabels(ax.get_xticklabels(), **xlab_kws)

    g.tight_layout()
    return g

def draw_heatmap(*args, **kwargs):
    data = kwargs.pop('data')
    yind = kwargs.pop('yind')
    d = data.pivot(index=args[1], columns=args[0], values=args[2])
    for y in yind:
        if y not in d.index:
            d.loc[y] = np.NaN
    d = d.sort_index()
    sns.heatmap(d, **kwargs)


def row_col_facet(data, func_args, func, row, col, row_order=None, col_order=None, **kwargs):
    if col_order: 
            cols = col_order
    else:
        cols = data[col].unique()
    if row_order:
        rows = row_order
    else:
        rows = data[row].unique()
    fg, axs = plt.subplots(len(rows), len(cols), figsize=(5*(len(cols) + 1), 5*len(rows)))
    for i in range(len(rows)):
        cax=None
        for j in range(len(cols)):
            plt_data = data[(data[row] == rows[i]) & (data[col] == cols[j])]
            cbar = j == (len(cols)- 1)
            if cbar:
                cax = axs[i][len(cols) - 1].inset_axes([1.05, 0, 0.05,1])
            func(*func_args, data=plt_data,ax=axs[i][j],
                                cbar=cbar,cbar_ax=cax,**kwargs)                            
            axs[i][j].set_title(f'{rows[i]}, {cols[j]}')
    return fg

def col_facet(data, func_args, func, col, col_wrap=None,col_order=None, **kwargs):
    if col_order:
        cols = col_order
    else:
        cols = data[col].unique()
    if col_wrap:
        fg, axs = plt.subplots((len(cols) // col_wrap) + 1,col_wrap, figsize=(5*col_wrap, 5*((len(cols) // col_wrap) + 1)))
    else:
        fg, axs = plt.subplots(1, len(cols), figsize=(5*len(cols), 5))
    axs = axs.flatten()
    cax=None
    for j in range(len(cols)):
        plt_data = data[(data[col] == cols[j])]
        cbar = j == (len(cols)-1)
        if cbar:
            cax = axs[len(cols) - 1].inset_axes([1.05, 0, 0.05,1])
        func(*func_args,data=plt_data,ax=axs[j],
                        cbar=cbar,cbar_ax=cax,**kwargs)                            
        axs[j].set_title(f'{cols[j]}')
    return fg

def row_facet(data, func_args, func, row,row_order=None, **kwargs):
    if row_order:
        rows = row_order
    else:
        rows = data[row].unique()    
    fg, axs = plt.subplots(len(rows), 1, figsize=(5, 5*len(rows)))
    for i in range(len(rows)):
        plt_data = data[(data[row] == rows[i])]
        cbar = True
        func(*func_args, data=plt_data,ax=axs[i],
                        cbar=cbar,**kwargs)                            
        axs[i].set_title(f'{rows[i]}')
    return fg

def single_facet(data, func_args, func,**kwargs):
    fg, ax = plt.subplots(1, 1, figsize=(5,5))
    cax = ax.inset_axes([1.05, 0, 0.05,1])
    func(*func_args, data=data,ax=ax,
                    cbar=True,cbar_ax=cax,**kwargs) 
    return fg

def draw_facet(data, func_args,  func, row=None, col=None, col_wrap=None, col_order=None, row_order=None,**kwargs):
    if row and col:
        cols = col_order if col_order else data[col].unique()
        rows = row_order if row_order else data[row].unique()
        yind = data[func_args[1]].unique()
        if len(rows) == 1 and len(cols) == 1:
            fg = single_facet(data, func_args, func,yind=yind, **kwargs)
            fg.axes[0][0].set_title(f'{rows[0]}, {cols[0]}')
        elif len(rows) == 1:
            fg = col_facet(data[data[row] == rows[0]], func_args, func, col, col_wrap=col_wrap, 
                           col_order=col_order,yind=yind,**kwargs)
        elif len(cols) == 1:
            fg = row_facet(data[data[col] == cols[0]], func_args, func, row, row_order,yind=yind, **kwargs)
        else:
            fg = row_col_facet(data, func_args, func, row, col,
                               row_order=row_order, col_order=col_order, yind=yind,**kwargs)
    elif row:
        rows = data[row].unique()
        if len(rows) == 1:
            fg = single_facet(data, func_args, func,yind=yind, **kwargs)
            fg.axes[0][0].set_title(f'{rows[0]}')
        else:
            fg = row_facet(data, func_args, func, row, row_order=row_order,yind=yind,**kwargs)
    elif col:
        cols = data[col].unique()
        if len(cols) == 1:
            fg = single_facet(data, func_args, func,yind=yind, **kwargs)
            fg.axes[0].set_title(f'{cols[0]}')
        else:
            fg = col_facet(data, func_args, func, col, col_wrap=col_wrap, 
                           col_order=col_order,**kwargs)
    else:
        fg = single_facet(data, func_args, func,yind=yind, **kwargs)
    fg.tight_layout()
    return fg

def plot_pair_strip(data, x=None, y=None, row=None, col=None, col_order=None, connect=False,
                   color=None, row_order=None,ylim=[0,10], xlim=[-0.5, 1.5], order=None, dodge=False, jitter=False,
                   hue=None, palette=None, hue_order=None, ylabel='', xlabel='',xlab_kws={}, aspect=0.5,legend=False,**kwargs):
    g = sns.catplot(data, y=y,x=x, order=order,kind="strip", hue=hue, 
                    palette=palette, hue_order=hue_order, col=col, row=row,
                    jitter=False, aspect=aspect, legend=legend,col_order=col_order, 
                    row_order=row_order, color=color, **kwargs)
    if g.axes_dict:
        for name, ax in g.axes_dict.items():
            if type(name) == tuple:
                cdfilt = data[(data[row] == name[0]) & (data[col] == name[1])]
                ax.set_title(f'{name[0]}, {name[1]}')
            else:
                if col:
                    cdfilt = data[data[col] == name]
                elif row:
                    cdfilt = data[data[row] == name]
            if connect:
                connect_lines(ax,cdfilt, hue=hue, order=order, x=x, y=y,hue_order=hue_order)                    
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.tick_params(labelbottom=True, labeltop=False)
            ax.set_xticklabels(ax.get_xticklabels(), **xlab_kws)
            ax.set_title(name)
    else:
        ax = g.axes[0][0]
        if connect:
                connect_lines(ax,cdfilt, hue=hue, order=order, x=x, hue_order=hue_order)                    
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.tick_params(labelbottom=True, labeltop=False)
        ax.set_xticklabels(ax.get_xticklabels(), **xlab_kws)
        name  = col if col else ''
        name  =  name + row if row else name
        ax.set_title(name)
    return g

def style_retrain(g, xlim, ylim, xticklab, xlab, ylab):
    g.axes[0][0].set_xticklabels(xticklab)
    xticks = g.axes[0][0].xaxis.get_major_ticks()
    if xlim[1] > 48:
        xticks[8].set_visible(False)
    if xlim[1] > 96:
        xticks[14].set_visible(False)
    xticks = g.axes[0][0].xaxis.get_minor_ticks()
    vis = []
    if xlim[1] > 48:
        vis = [13, 14, 15, 16]
        g.axes[0][0].text(60, ylim[0]-0.06*(ylim[1]-ylim[0]), '//', ha='center', va='bottom', fontsize=24)
        g.axes[0][0].text(96, ylim[0]-.2*(ylim[1]-ylim[0]), '5-8d test', ha='center')
    if xlim[1] > 96:
        vis = vis + [25, 26, 27, 28]
        g.axes[0][0].text(144-12, ylim[0]-0.06*(ylim[1]-ylim[0]), '//', ha='center', va='bottom', fontsize=24)
        g.axes[0][0].text(168, ylim[0]-.2*(ylim[1]-ylim[0]), '21-24d test', ha='center')
    for i in vis:
            xticks[i].set_visible(False)
    g.axes[0][0].text(12, ylim[0]-.2*(ylim[1]-ylim[0]), 'Train', ha='center')
    g.axes[0][0].set_xlabel(xlab, labelpad=28)
    g.axes[0][0].set_ylabel(ylab)


def plot_retrain(c, means, perf, tot, xlim, xticklab,xlab='Time (hr)', lckylab='Lick Frequency (Hz)', 
                 perfylab="Performance\n$L_s - L_b$", totylab='Number of trials',
                 lckylim=[-0.1, 10], perfylim=[-6,6], totylim=[0, 400], aspect=2, **kwargs):
    g1 = None
    if not means.empty:
        g1 = plot_ant_lickfreq(means[(means['cond'] == c) & (means['test_time'] == 'train')],
                            y='lick',aspect=aspect,lw=1, ms=7,col='cond', title=c, 
                            legend=False, xlim=xlim, ylim=lckylim, errorbar='se', 
                            err_style='bars', err_kws={"lw":0.5},wspace=0.15,style=False, **kwargs)
        plt_data = means[(means['cond'] == c) & (means['test_time'] == '5d test')]
        if not plt_data.empty:
            sns.lineplot(plt_data,
                        y='lick',lw=1, ms=7,ax=g1.axes[0][0],x='Time (hr)',
                        errorbar='se', err_style='bars', err_kws={"lw":0.5},hue='stimulus', 
                        hue_order=['blank', 'stimulus'], palette=['red', 'green'], 
                        mec=None, marker='o', legend=False)
        plt_data = means[(means['cond'] == c) & (means['test_time'] == '21d test')]
        if not plt_data.empty:
            sns.lineplot(plt_data,
                        y='lick',lw=1, ms=7,ax=g1.axes[0][0],x='Time (hr)',
                        errorbar='se', err_style='bars', err_kws={"lw":0.5},hue='stimulus', 
                        hue_order=['blank', 'stimulus'], palette=['red', 'green'], 
                        mec=None, marker='o',legend=False)
        style_retrain(g1, xlim, lckylim, xticklab, xlab, lckylab)
    g2 = None
    if not perf.empty:
        g2 = plot_ant_perf(perf[(perf['cond'] == c) & (perf['test_time'] == 'train')],
                        ylabel=perfylab, aspect=aspect,lw=1, ms=7,
                        title='',xlim=xlim,col='condition',legend=True, ylim=perfylim, 
                        errorbar='se', err_style='bars', err_kws={"lw":0.5},
                        wspace=0.25,**kwargs)
        sns.lineplot(perf[(perf['cond'] == c) & (perf['test_time'] == '5d test')],
                    y='lick',lw=1, ms=7,ax=g2.axes[0][0],x='Time (hr)', errorbar='se', 
                    err_style='bars', err_kws={"lw":0.5}, mec=None, marker='o', 
                    color='k', legend=False)
        sns.lineplot(perf[(perf['cond'] == c) & (perf['test_time'] == '21d test')],
                    y='lick',lw=1, ms=7,ax=g2.axes[0][0],x='Time (hr)', errorbar='se', 
                    err_style='bars', err_kws={"lw":0.5},color='k', mec=None, 
                    marker='o', legend=False)
        style_retrain(g2, xlim, perfylim, xticklab,xlab, perfylab)
    g3 = None
    if not tot.empty:
        g3 = plot_trial_hr(tot[(tot['cond'] == c) & (tot['test_time'] == 'train')],
                        xlim=xlim, col='condition',fill=False,ylim=totylim, 
                        aspect=aspect*1.5, title='', ms=7, lw=1,err_kws={"lw":0.5}, 
                        ylabel=totylab,wspace=0.04,**kwargs)
        sns.lineplot(tot[(tot['cond'] == c) & (tot['test_time'] == '5d test')],
                    y='trial no',lw=1, ms=7,ax=g3.axes[0][0],x='Time (hr)', 
                    errorbar='se', err_style='bars', err_kws={"lw":0.5},color='k', 
                    mec=None, marker='o', legend=False)
        sns.lineplot(tot[(tot['cond'] == c) & (tot['test_time'] == '21d test')],
                    y='trial no',lw=1, ms=7,ax=g3.axes[0][0],x='Time (hr)',
                    errorbar='se', err_style='bars', err_kws={"lw":0.5},color='k', 
                    mec=None, marker='o', legend=False)
        style_retrain(g3, xlim, totylim, xticklab, xlab, totylab)
        for i in range(len(g3.axes[0][0].lines)):
            if i % 2 == 0:
                    g3.axes[0][0].fill_between(g3.axes[0][0].lines[i].get_data()[0], g3.axes[0][0].lines[i].get_data()[1], 
                                            color='k', alpha=0.3)
    return (g1, g2, g3)