import seaborn as sns #v0.13.0
import matplotlib as mpl #v3.8.1
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.dates as mdates

#import matplotlib.spines as spines
from matplotlib.patches import Rectangle



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
             suptitle=None, title=None, ymajmult=2, yminmult=1):
    style_axes(g, xlabel, ylabel, xlim, ylim, 12, 4, ymajmult, yminmult, 
               hspace=hspace, wspace=wspace, suptitle=suptitle, title=title)
    if g.axes_dict:
        for (name,ax) in g.axes_dict.items():
            ax.add_patch(Rectangle((0, ylim[0]), xlim[1], ylim[1] - ylim[0], color="#CEE1F3", alpha=0.8, 
                                zorder=0, fill=True))        
    else:
        ax = g.figure.axes[0]
        ax.add_patch(Rectangle((0, ylim[0]), xlim[1], ylim[1] - ylim[0], color="#CEE1F3", alpha=0.8, 
                                zorder=0, fill=True))     

def style_ms(g, xlabel, ylabel, xlim, ylim, hspace=0.45, wspace=0.25, 
             suptitle=None):
    style_axes(g, xlabel, ylabel, xlim, ylim, 500, 200, 2, 1, hspace=hspace,
               wspace=wspace, suptitle=suptitle)
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
                    err_kws=err_kws, sharex=sharex, sharey=sharey,
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
                  xlabel="Time (ms)", hspace=0.45, wspace=0.25, suptitle=None, **kwargs):
    
    g = lineplot(data, x=x, y=y, hue=hue, col_wrap=col_wrap, col=col,
                 palette=palette, color=color, hue_order=hue_order, 
                 legend=legend, **kwargs)
    style_ms(g, xlabel, ylabel, xlim, ylim, hspace=hspace, wspace=wspace, 
             suptitle=suptitle)
    
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

def plot_trial_bar(data, x="condition", y="trial no", aspect=1, col="Day", 
                   color="#3B3838", palette=None, col_wrap=None, hue_order=None,
                   hue=None, xlim=[-0.5, 2.5], ylim=[0, 600], xlabel=None, hspace=0.45, wspace=0.25,
                   ylabel="Number of Trials", width=0.4, ymajmult=100, yminmult=50,**kwargs):
    g = barplot(data, x=x, y=y, aspect=aspect, col=col, color=color, 
                palette=palette, hue_order=hue_order, hue=hue, 
                col_wrap=col_wrap, width=width, **kwargs)
    style_axes(g, xlabel, ylabel, xlim, ylim, None, None, ymajmult=ymajmult, 
               yminmult=yminmult, hspace=hspace, wspace=wspace)

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
    for name, ax in g.axes_dict.items():
        ax.spines[['right']].set_visible(True)
        ax.spines[['left']].set_visible(False)
        ax.xaxis.set_visible(False)
        ax.yaxis.tick_right()

    if fill: style_trial(g)
    return g