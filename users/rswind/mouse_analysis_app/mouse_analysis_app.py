# -*- coding: utf-8 -*-
"""
Created 2024-05-29
Last Edited 2024-07-02

@author: Rachel Swindell
"""

import os
from pathlib import Path
import re
import datetime
import math
from collections import Counter

import tkinter as tk
from tkinter import ttk, filedialog, StringVar, BooleanVar

import pandas as pd
import numpy as np

import matplotlib as mpl
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle
from matplotlib.collections import LineCollection

import seaborn as sns #v0.13.0

class WrappingLabel(tk.Label):
    def __init__(self, master=None, **kwargs):
        tk.Label.__init__(self, master, **kwargs)
        self.bind('<Configure>', lambda e: self.config(wraplength=self.winfo_width()))

class ResizeEqualLabelFrame(ttk.LabelFrame):
    def __init__(self, master=None,cols=0, rows=0, **kwargs):
        ttk.LabelFrame.__init__(self, master, **kwargs)
        for i in range(cols):
            self.grid_columnconfigure(i, weight=1)
        for i in range(rows):
            self.grid_rowconfigure(i, weight=1)

class EnterButton(ttk.Button):
    def __init__(self, master=None,**kwargs):
        ttk.Button.__init__(self, master, **kwargs)
        self.bind('<Return>', lambda e: self.invoke())
        
class EnterCheckbutton(ttk.Checkbutton):
    def __init__(self, master=None,**kwargs):
        ttk.Checkbutton.__init__(self, master, **kwargs)
        self.bind('<Return>', lambda e: self.invoke())

class CustomToolbar(NavigationToolbar2Tk):
    def __init__(self, canvas_, parent_, **kwargs):
        self.toolitems = (
        ('Home', 'Reset original view', 'home', 'home'),
        ('Back', 'Back to  previous view', 'back', 'back'),
        ('Forward', 'Forward to next view', 'forward', 'forward'),
        (None, None, None, None),
        ('Pan', 'Pan axes with left mouse, zoom with right', 'move', 'pan'),
        ('Zoom', 'Zoom to rectangle', 'zoom_to_rect', 'zoom'),
        (None, None, None, None),
        ('Subplots', 'Configure subplots', 'subplots', 'configure_subplots'),
        ('Advanced', 'Advanced figure control', 'advanced', 'advanced_window'),
        (None, None, None, None),
        ('Save', 'Save the figure', 'filesave', 'save_figure'),
        )
        NavigationToolbar2Tk.__init__(self, canvas_, parent_, **kwargs)

    def advanced_window(self):
        AdvancedWindow(self)

class MplEmbedPlot():
    def __init__(self,fig=None, master=None,**kwargs):
        self.frame = ttk.Frame(master, **kwargs)
        self.canvas = FigureCanvasTkAgg(fig, master=self.frame)
        self.canvas.draw()
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.frame, pack_toolbar=False)
        #self.toolbar = CustomToolbar(self.canvas, self.frame, pack_toolbar=False)
        self.toolbar.update()
        self.fontlabel = tk.Label(text='Font size')
        #self.fontentry = ttk.Entry()

        #self.fontentry.bind('<FocusOut>', self.update_font)
        #self.fontentry.bind('<Return>', self.update_font)
        self.frame.bind_all('<Escape>', lambda e: e.widget.focus_set())

        self.toolbar.pack(side=tk.BOTTOM, fill=tk.X)
        #self.fontlabel.pack(side=tk.BOTTOM)
        #self.fontentry.pack(side=tk.BOTTOM, fill=tk.X)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        

    def update_font(self, e):
        size = int(self.fontentry.get())
        for ax in self.canvas.figure.axes:
            ax.tick_params(axis='both', which='major', labelsize=size)
            ax.xaxis.label.set_size(size)
            ax.yaxis.label.set_size(size)
            ax.title.label.set_size(size)
        self.canvas.figure.suptitle.set_size(size)
        
class AdvancedWindow(tk.Toplevel):
    def __init__(self, master=None,**kwargs):
        tk.Toplevel.__init__(self, master, **kwargs)
        self.frame = ResizeEqualLabelFrame(self)
        self.title('Advanced Plot Adjustment')
        self.validate_number_wrapper = (self.frame.register(self.validate_number), '%P')
        # row
        self.rllab = ttk.Label(self.frame,text="Row names")
        self.rlval = StringVar()
        self.rlen = ttk.Entry(self.frame, textvariable=self.rlval)

        # col
        self.clab = ttk.Label(self.frame,text="Number of columns")
        self.cval = StringVar()
        self.cent = ttk.Entry(self.frame,textvariable=self.cval)
        self.cllab = ttk.Label(self.frame,text="Column names")
        self.clvar = StringVar()
        self.clen = ttk.Entry(self.frame, textvariable=self.clvar)

        # time start and end
        self.thslab = ttk.Label(self.frame,text="Time (hr) start (for instantaneous plots only)")
        self.thelab = ttk.Label(self.frame,text="Time (hr) end (for instantaneous plots only)")
        self.thsvar = StringVar()
        self.thevar = StringVar()
        self.thsen = ttk.Entry(self.frame, textvariable=self.thsvar, validatecommand=self.validate_number_wrapper)
        self.theen = ttk.Entry(self.frame, textvariable=self.thevar, validatecommand=self.validate_number_wrapper)
        
        # x and y lim
        self.xminlab = ttk.Label(self.frame,text="X min")
        self.xmaxlab = ttk.Label(self.frame,text="X max")
        self.xminvar = StringVar()
        self.xmaxvar = StringVar()
        self.xminen = ttk.Entry(self.frame, textvariable=self.xminvar, 
                              validatecommand=self.validate_number_wrapper, 
                              invalidcommand=(lambda: self.xminvar.set('')))
        self.xmaxen = ttk.Entry(self.frame, textvariable=self.xmaxvar, 
                              validatecommand=self.validate_number_wrapper,
                              invalidcommand=(lambda: self.xmaxvar.set('')))
        
        self.yminlab = ttk.Label(self.frame,text="Y min")
        self.ymaxlab = ttk.Label(self.frame,text="Y max")
        self.yminvar = StringVar()
        self.ymaxvar = StringVar()
        self.yminen = ttk.Entry(self.frame, textvariable=self.yminvar, 
                              validatecommand=self.validate_number_wrapper, 
                              invalidcommand=(lambda: self.yminvar.set('')))
        self.ymaxen = ttk.Entry(self.frame, textvariable=self.ymaxvar, 
                              validatecommand=self.validate_number_wrapper,
                              invalidcommand=(lambda: self.ymaxvar.set('')))

        # Which Days to plot
        self.dayslab = ttk.Label(self.frame,text='Days (for nth part only)')
        self.daysvar = StringVar()
        self.daysen = ttk.Entry(self.frame, textvariable=self.daysvar)
        
        # view changes
        # font size
        self.fontsizelab = ttk.Label(self.frame,text='fontsize')
        self.fsvar = StringVar()
        self.fsen = ttk.Entry(self.frame,textvariable=self.fsvar, 
                              validatecommand=self.validate_number_wrapper,
                              invalidcommand=(lambda: self.fsvar.set('')))

        # number of ticks
        self.numtkxlab = ttk.Label(self.frame, text="Number of x ticks")
        self.nxtvar = StringVar()
        self.numbtkxen = ttk.Entry(self.frame,textvariable=self.nxtvar, 
                              validatecommand=self.validate_number_wrapper,
                              invalidcommand=(lambda: self.nxtvar.set('')))
        self.numtkylab = ttk.Label(self.frame, text='Number of y ticks')
        self.nytvar = StringVar()
        self.numbtkyen = ttk.Entry(self.frame, textvariable=self.nytvar, 
                              validatecommand=self.validate_number_wrapper,
                              invalidcommand=(lambda: self.nytvar.set('')))
        
        # labels and titles
        self.titlab = ttk.Label(self.frame, text='Title')
        self.titvar = StringVar()
        self.titen = ttk.Entry(self.frame, textvariable=self.titvar)
        self.xlablab = ttk.Label(self.frame, text='xlabel')
        self.xlabvar = StringVar()
        self.xlabent = ttk.Entry(self.frame, textvariable=self.xlabvar)
        self.xbot = BooleanVar(value=False)
        self.xbotonly = EnterCheckbutton(self.frame, text='xlabel on bottom plots only', 
                                         onvalue=True, offvalue=False, variable=self.xbot)

        self.ylablab = ttk.Label(self.frame, text='ylabel')
        self.ylabvar = StringVar()
        self.ylabent = ttk.Entry(self.frame, textvariable=self.ylabvar)
        self.ybot = BooleanVar(value=False)
        self.ybotonly = ttk.Checkbutton(self.frame, text='ylabel on left plots only', 
                                        onvalue=True, offvalue=False, variable=self.ybot)

        self.figwidlab = ttk.Label(self.frame, text='Figure width (inches)')
        self.fwvar = StringVar()
        self.fwen = ttk.Entry(self.frame, textvariable=self.fwvar, 
                              validatecommand=self.validate_number_wrapper,
                              invalidcommand=(lambda: self.fwvar.set('')))

        self.figlenlab = ttk.Label(self.frame, text='Figure length (inches)')
        self.flvar = StringVar()
        self.flen = ttk.Entry(self.frame, textvariable=self.flvar, 
                              validatecommand=self.validate_number_wrapper,
                              invalidcommand=(lambda: self.flvar.set('')))
        
        # line features
        # lw, ms, mec,
        # change ms to s for pair_plot_strip

        self.plotbutt = EnterButton(self.frame,text='Plot', command=self.generate_plots)

        self.frame.grid(row=0, column=0, sticky='nsew')
        self.rllab.grid(row=1, column=0, padx=2, pady=2, sticky='nsew')
        self.rlen.grid(row=1, column=1, padx=2, pady=2, sticky='nsew')

        self.clab.grid(row=2, column=0, padx=2, pady=2, sticky='nsew')
        self.cent.grid(row=2, column=1, padx=2, pady=2, sticky='nsew')
        self.cllab.grid(row=3, column=0, padx=2, pady=2, sticky='nsew')
        self.clen.grid(row=3, column=1, padx=2, pady=2, sticky='nsew')

        self.thslab.grid(row=4, column=0, padx=2, pady=2, sticky='nsew')
        self.thsen.grid(row=4, column=1, padx=2, pady=2, sticky='nsew')
        self.thelab.grid(row=5, column=0, padx=2, pady=2, sticky='nsew')
        self.theen.grid(row=5, column=1, padx=2, pady=2, sticky='nsew')

        self.xminlab.grid(row=6, column=0, padx=2, pady=2, sticky='nsew')
        self.xminen.grid(row=6, column=1, padx=2, pady=2, sticky='nsew')        
        self.xmaxlab.grid(row=7, column=0, padx=2, pady=2, sticky='nsew')
        self.xmaxen.grid(row=7, column=1, padx=2, pady=2, sticky='nsew')

        self.yminlab.grid(row=8, column=0, padx=2, pady=2, sticky='nsew')
        self.yminen.grid(row=8, column=1, padx=2, pady=2, sticky='nsew')        
        self.ymaxlab.grid(row=9, column=0, padx=2, pady=2, sticky='nsew')
        self.ymaxen.grid(row=9, column=1, padx=2, pady=2, sticky='nsew')

        self.dayslab.grid(row=10, column=0, padx=2, pady=2, sticky='nsew')
        self.daysen.grid(row=10, column=1, padx=2, pady=2, sticky='nsew')

        self.fontsizelab.grid(row=11, column=0, padx=2, pady=2, sticky='nsew')
        self.fsen.grid(row=11, column=1, padx=2, pady=2, sticky='nsew')

        self.numtkxlab.grid(row=12, column=0, padx=2, pady=2, sticky='nsew')
        self.numbtkxen.grid(row=12, column=1, padx=2, pady=2, sticky='nsew')
        self.numtkylab.grid(row=13, column=0, padx=2, pady=2, sticky='nsew')
        self.numbtkyen.grid(row=13, column=1, padx=2, pady=2, sticky='nsew')

        self.titlab.grid(row=14, column=0, padx=2, pady=2, sticky='nsew')
        self.titen.grid(row=14, column=1, padx=2, pady=2, sticky='nsew')

        self.xlablab.grid(row=15, column=0, padx=2, pady=2, sticky='nsew')
        self.xlabent.grid(row=15, column=1, padx=2, pady=2, sticky='nsew')
        self.xbotonly.grid(row=16, column=0, padx=2, pady=2, sticky='nsew')
        
        self.ylablab.grid(row=17, column=0, padx=2, pady=2, sticky='nsew')
        self.ylabent.grid(row=17, column=1, padx=2, pady=2, sticky='nsew')
        self.ybotonly.grid(row=18, column=0, padx=2, pady=2, sticky='nsew')

        self.figwidlab.grid(row=19, column=0, padx=2, pady=2, sticky='nsew')
        self.fwen.grid(row=19, column=1, padx=2, pady=2, sticky='nsew')

        self.figlenlab.grid(row=20, column=0, padx=2, pady=2, sticky='nsew')
        self.flen.grid(row=20, column=1, padx=2, pady=2, sticky='nsew')
        
        self.plotbutt.grid(row=21, column=0, padx=2, pady=2, sticky='nsew')
        self.validate_number_wrapper = (self.frame.register(self.validate_number), '%P')

    def validate_number(self, newval):
        valid = re.match('^[0-9]*\.?[0-9]*$', newval) is not None
        return valid
    
    def generate_ant_lf(self, **kwargs):
        ant_means = self.data['fixed means filtered']
        ant_means['Time (hr)'] = ant_means['Time (hr)'] + 2
        g = self.master.plot_ant_lickfreq(ant_means, aspect=2, ms=7, **kwargs)
        g.figure.subplots_adjust(left=0.05, bottom=0.05, top=0.95, right=0.95)
        self.master.lfcan.canvas = FigureCanvasTkAgg(g, master=self.master.lfcan.frame)
        return g
            
    def generate_ant_perf(self, **kwargs):
        w = kwargs.pop('w', '')
        h = kwargs.pop('h', '')
        perf = self.data['fixed performance filtered']    
        perf['Time (hr)'] = perf['Time (hr)'] + 2
        g = self.master.plot_ant_perf(perf, aspect=2,ms=7, **kwargs)
        g.figure.subplots_adjust(left=0.05, bottom=0.05, top=0.95, right=0.95)
        if w != '' and h != '':
            g.figure.set_size_inches(w, h)
        self.master.perfcan.canvas = FigureCanvasTkAgg(g, master=self.master.perfcan.frame)

    def genereate_fullckhr(self, **kwargs):
        fontsize = int(kwargs.pop("fontsize", '8'))
        w = kwargs.pop("w", '10') 
        h = kwargs.pop("h", '')
        hrs = kwargs.pop('start_hour', '')
        hre = kwargs.pop('end_hour', '')
        means_full = self.data['inst means filtered']
        if hrs == '': hrs = means_full['Time (hr)'].min()
        else: hrs = int(hrs)
        if hre == '': hre = means_full['Time (hr)'].min()
        else: hre = int(hre)
        for c in means_full['condition'].unique():
            plt_data = means_full[(means_full['condition'] == c) & (means_full['Time (hr)'] >= hrs)
                                   & (means_full['Time (hr)'] <= hre)]
            g = self.master.plot_lickfreq(plt_data,col_wrap=6, ms=1, lw=0.7, col='Time (hr)',
                            wspace=0.2, hspace=0.2, **kwargs)
            g.tick_params(labelsize=8, pad=0.5)
            i=0
            for n, ax in g.axes_dict.items():
                ax.set_xlabel(ax.get_xlabel(), fontsize=fontsize)
                ax.set_ylabel(ax.get_ylabel(), fontsize=fontsize)
                ax.set_title(f'Time (hr):{ax.get_title()}', fontsize=fontsize)
                if i % 6 != 0:
                    ax.set_ylabel('')
                if i < len(g.axes)-6:
                    ax.set_xlabel('')
                ax.tick_params(axis='x', which='major', labelrotation=45)
                i +=1
            nrow = math.ceil(len(g.axes)/6)
            h1 = nrow*1.5 if h == '' else int(h)
            g.figure.set_size_inches(w, h1)
            g.legend.set_in_layout(True)
            g.figure.tight_layout()

            g.legend.set_title("Stimulus Type", prop=mpl.font_manager.FontProperties(size=fontsize))
            for t in g.legend.get_texts():
                t.set_size(fontsize)
            g.legend.set(loc='center right',bbox_to_anchor=(1.01, 0.5))
            g.figure.subplots_adjust(left=0.1, bottom=0.15,right=0.9)
            self.master.fllfplots[c].canvas = FigureCanvasTkAgg(g, master=self.master.perfcan.frame)

    def generate_fullpfhr(self, **kwargs):
        fontsize = int(kwargs.pop("fontsize", '8'))
        w = kwargs.pop("w", '10') 
        h = kwargs.pop("h", '')
        perf = self.data['inst performance filtered']
        for c in perf['condition'].unique():
            g = self.master.plot_perf(perf[perf['condition'] == c],  
                        col="Time (hr)", col_wrap=6, hue=None, color='k',suptitle=c,
                        ms=2, lw=0.7, wspace=0.2,hspace=0.2)
            g.tick_params(labelsize=8, pad=0.5)
            i=0
            nrow = math.ceil(len(g.axes)/6)
            for n, ax in g.axes_dict.items():
                ax.set_xlabel(ax.get_xlabel(), fontsize=8)
                ax.set_ylabel(ax.get_ylabel(), fontsize=8)
                ax.set_title(f'Time (hr):{ax.get_title()}', fontsize=8)
                if i % 6 != 0:
                    ax.set_ylabel('')
                if i < len(g.axes)-6:
                    ax.set_xlabel('')
                ax.tick_params(which='major', axis='x', labelrotation=45)
                i +=1
            g.figure.set_size_inches(6*1.5, nrow*1.5)
            g.tight_layout()
            g.figure.subplots_adjust(left=0.05, bottom=0.15)

    def generate_htmplf(self):
        if self.htmplf:
            try:
                self.data['inst means filtered']
            except KeyError:
                errdisp = tk.Toplevel()
                errlab = tk.Label(errdisp, text='No instantaneous analysis')
                errlab.grid()
                return
            gpcols = ['condition', 'Time (hr)', 'Time (ms)','stimulus']
            means_full_means = self.get_means(self.data['inst means filtered'], gpcols)
            htmp_data = means_full_means.sort_values(['Time (hr)', 'Time (ms)'])
            htmp_data = htmp_data[(htmp_data['Time (ms)'] < 2000) & (htmp_data['Time (ms)'] > -700)]
            htmp_args=['Time (ms)', 'Time (hr)', 'lick']
            fg = self.draw_facet(htmp_data, htmp_args, self.draw_heatmap, col='stimulus', row='condition',
                    col_order=['stimulus', 'blank'],square=True,
                    cmap='viridis', vmin=0, vmax=10, xticklabels=3)
            fg.suptitle('Instantaneous licking', y=0.98)
            fg.subplots_adjust(left=0.05, bottom=0.15, top=0.9, right=0.95, wspace=0, hspace=0.75)
            if self.disp:
                self.fdisp = MplEmbedPlot(fg.figure, master=self.plt_container)
                self.plt_container.add(self.fdisp.frame, text='Inst Lick Freq Heatmap')
            if self.save:
                fg.savefig(self.savedir+"\\lckfreqheatmap.png", transparent=True)

    def generate_htmppf(self):
        if self.htmppf:
            try:
                self.data['inst performance filtered']
            except KeyError:
                errdisp = tk.Toplevel()
                errlab = tk.Label(errdisp, text='No instantaneous analysis')
                errlab.grid()
                return
            gpcols = ['condition', 'Time (hr)', 'Time (ms)']
            perf_full_means = self.get_means(self.data['inst performance filtered'], gpcols)
            htmp_data = perf_full_means.sort_values(['Time (hr)', 'Time (ms)'])
            htmp_data = htmp_data[(htmp_data['Time (ms)'] < 2000) & (htmp_data['Time (ms)'] > -700)]
            htmp_args=['Time (ms)', 'Time (hr)', 'lick']
            fg = self.draw_facet(htmp_data, htmp_args, self.draw_heatmap, col='condition',square=True,
                    cmap='coolwarm', vmin=-6, vmax=6, col_wrap=6, xticklabels=3)
            fg.suptitle('Instantaneous performance', y=0.98)
            fg.subplots_adjust(left=0.05, bottom=0.1, top=0.9, right=0.95, wspace=0.25, hspace=0.25)
            if self.disp:
                self.pddisp = MplEmbedPlot(fg.figure, master=self.plt_container)
                self.plt_container.add(self.pdisp.frame, text='Inst Perf Heatmap')
            if self.save:
                fg.savefig(self.savedir+"\\perfheatmap.png", transparent=True)
 
    def generate_nthpart(self):
        if self.nthpart:
            try:
                last20 = self.data['part means']
            except KeyError:
                errdisp = tk.Toplevel()
                errlab = tk.Label(errdisp, text='No part analysis')
                errlab.grid()
                return
            last20["Day"] = last20["Day"].apply(self.day_to_label)
            g = self.plot_pair_strip(last20, y="lick",x="stimulus", order=["stimulus", "blank"],
                            hue="stimulus", palette=["green", "red"], ylabel='Lick frequency (Hz)',
                            hue_order=["stimulus", "blank"],
                            row="condition", col="Day",jitter=False, s=7, aspect=0.5, legend=False, connect=True)
            g.figure.subplots_adjust(left=0.05, bottom=0.05, top=0.95, right=0.95)
            for (r, c), ax in g.axes_dict.items():
                ax.set_title(f'{r} {c}', fontsize=10)
            if self.disp:
                self.nthdisp = MplEmbedPlot(g.figure, master=self.plt_container)
                self.plt_container.add(self.nthdisp.frame, text='Nth Part')
            if self.save:
                g.savefig(self.savedir+"\\nthpart.png", transparent=True)

    def generate_trialshr(self):
        if self.trialshour:
            try:
                trials = self.data['fixed trial counts']
            except KeyError:
                errdisp = tk.Toplevel()
                errlab = tk.Label(errdisp, text='No fixed analysis')
                errlab.grid()
                return
            tots = self.sum_trials(trials, ['condition', 'animal', 'Time (hr)'],[], "trial no")
            xlim = [tots['Time (hr)'].min() - 2, tots['Time (hr)'].max() + 2]
            ymax = math.ceil(tots['trial no'].max()/100.0)*100
            g = self.plot_trial_hr(tots, xlim=xlim,ylim=[0, ymax], title=None, ms=7, lw=1,err_kws={"lw":0.5}, ylabel='Trial no')
            g.figure.subplots_adjust(left=0.05, bottom=0.05, top=0.95, right=0.95)
            if self.disp:
                self.trialsdisp = MplEmbedPlot(g.figure, master=self.plt_container)
                self.plt_container.add(self.trialsdisp.frame, text='Trials by hour')
            if self.save:
                g.savefig(self.savedir +"\\trialsbyhour.png", transparent=True)
    
    def day_to_label(self,day):
        day = int(day)
        if day < 0:
            return "ACC " + str(day + 1)
        return "SAT " + str(day + 1)
    
    def sum_trials(self,data, group, keep, key):
        '''Get total number of trials grouped by given group. Only use on aggregated
        trial data.
        
        Parameters:
        data --- data frame containing trial data. Must contain all columns in
                `group`, `keep`, and `key`.
        group --- list of columns to group by
        keep --- list of columns to keep
        key --- column containing trial counts'''
        tot = data.groupby(group)[key].sum().rename('trial no')
        keep = data.groupby(group).first()[keep]
        res = pd.concat([tot, keep], axis=1).reset_index()
        return res
   
    
    def generate_plots(self):
        curr_plot_name = self.master.master.plot_container.select()

        row  = None if self.rlval.get() == '' else self.rlval.get()
        numcols = None if self.cval.get() == '' else int(self.cval.get())
        col = None if self.clvar.get() == '' else self.clvar.get() 

        ths = self.thsvar.get()
        the = self.thevar.get()

        xmin = None if self.xminvar.get() == '' else int(self.xminvar.get())
        xmax = None if self.xmaxvar.get() == '' else int(self.xmaxvar.get())
        ymin = None if self.xminvar.get() == '' else int(self.xminvar.get())
        ymax = None if self.xmaxvar.get() == '' else int(self.xmaxvar.get())

        days = self.daysvar.get()

        fontsize = None if self.fsvar.get() == '' else self.fsvar.get()
        xticks = None if self.nxtvar.get() == '' else int(self.nxtvar.get())
        yticks = None if self.nytvar.get() == '' else int(self.nytvar.get())

        title = None if self.titvar.get() == '' else self.titvar.get()
        xlabel = None if self.xlabvar.get() == '' else self.xlabvar.get()
        xlabbot = self.xbot.get()
        ylablef = self.ybot.get()
        ylabel = None if self.ylabvar.get() == '' else self.ylabvar.get()
        w = None if self.fwvar.get() == '' else int(self.fwvar.get())
        h = None if self.flvar.get() == '' else int(self.flvar.get())
        match curr_plot_name:
            case "Fixed Lick Frequency":
                curr_plot = self.master.lfcan.figure
                old_ymin, old_ymax = curr_plot.get_ylim()
                old_xmin, old_xmax = curr_plot.get_xlim()
                ymin = old_ymin if ymin == '' else int(ymin)
                ymax = old_ymax if ymax == '' else int(ymax)
                xmin = old_xmin if xmin == '' else int(xmin)
                xmax = old_xmax if xmax == '' else int(xmax)

                g = self.master.generate_ant_lf(row=row, col=col, col_wrap=numcols, ylim=[ymin, ymax], xlim=[xmin,xmax], ylabel=ylabel, xlabel=xlabel,title=title )
                if w != None and h != None:
                    g.figure.set_size_inches(w, h)
                if g.axes_dict:
                    for ax in g.axes_dict.items:
                        if yticks:
                            ax.yaxis.set_major_locator(ticker.LinearLocator(yticks))
                            ax.yaxis.set_minor_locator(ticker.LinearLocator(4*(yticks-1)+1))
                        if xticks:
                            ax.xaxis.set_major_locator(ticker.LinearLocator(xticks))
                            ax.xaxis.set_minor_locator(ticker.LinearLocator(4*(xticks-1)+1))
                else:
                    ax = g.axes[0]
                    if yticks:
                        ax.yaxis.set_major_locator(ticker.LinearLocator(yticks))
                        ax.yaxis.set_minor_locator(ticker.LinearLocator(4*(yticks-1)+1))
                    if xticks:
                        ax.xaxis.set_major_locator(ticker.LinearLocator(xticks))
                        ax.xaxis.set_minor_locator(ticker.LinearLocator(4*(xticks-1)+1))
                    
                #self.master.lfcan.canvas.draw()
                return
            case "Fixed Lick Performance":
                curr_plot = self.master.perfcan.figure
                self.master.generate_ant_perf(row=row, col=col, col_wrap=numcols, fontsize=fontsize, title=title, xlabel=xlabel, ylabel=ylabel)
                return
            case "Instentaneous Performance":
                curr_plot = self.master.flperfcan
                self.master.genereate_fullperf()
                return
            case "Instentaneous Lick Freq":
                curr_plot = self.master.fllfplots
                self.master.genereate_fullckhr()
                return
            case "Inst Lick Freq Heatmap":
                curr_plot = self.master.fdisp.figure
                self.master.generate_htmplf()
                return
            case "Inst Perf Heatmap":
                curr_plot = self.master.htmppfdisp.figure
                self.master.generate_htmppf()
                return
            case "Nth Part":
                curr_plot = self.master.nthdisp.figure
                self.master.generate_nthpart()
                return
            case "Trials by hour":
                curr_plot = self.master.trialsdisp.figure
                self.master.generate_trialshr()
                return

        return

        #self.master.generate_plots()

class PlotWindow(tk.Toplevel):
    def __init__(self, master=None, **kwargs):
        tk.Toplevel.__init__(self, **kwargs)
        self.geometry('800x600+100+100')
        self.plt_container = ttk.Notebook(self)
        self.plt_container.pack()

        self.antlf = master.ant_lf_var.get()
        self.antperf = master.ant_perf_var.get()
        self.instlf = master.fulllckhr.get()
        self.instperf = master.fullpfhr.get()
        self.htmplf = master.htmplck.get()
        self.htmppf = master.htmppf.get()
        self.trialshour = master.trialsbyhour.get()
        self.nthpart = master.last20val.get()

        self.loaddir = master.loadfilevar.get()
        self.savedir = master.savefilevar.get()
        self.save = master.savevar.get()
        self.disp = master.dispvar.get()
        
        self.dtype_map = {'condition':str,'animal':str, 'stimulus':str, 'lick':float, 'Strain':str, 'Cage':str, 'type':str, 'Time (ms)':float, 'Time (hr)':float, 'Day':float, 'water':str}

        # load data
        data = dict()
        for type_name in os.listdir(self.loaddir):
            if 'data' not in type_name and os.path.isfile(self.loaddir + '\\' + type_name):
                splt_name = type_name.split('.')[0].split('_')
                splt_name = [re.search('[a-z]+',splt_name[0])[0]] + splt_name[2:]
                data[' '.join(splt_name)] = pd.read_csv(self.loaddir + '\\' + type_name, dtype=self.dtype_map, header=0)
        self.data = data

        self.generate_ant_perf()
        self.generate_ant_lf()
        self.generate_fullpfhr()
        self.genereate_fullckhr()
        self.generate_htmplf()
        self.generate_htmppf()
        self.generate_nthpart()
        self.generate_trialshr()

    ##########################
    # Functions for Plotting #
    ##########################

    def get_means(self, df, gpcols):
        means = df.groupby(gpcols)['lick'].mean()
        kpcols = df.columns.drop(gpcols + ['lick'])
        kp  = df.groupby(gpcols)[kpcols].first()
        means = pd.concat([means, kp], axis=1).reset_index()
        return means

    def style_axes_helper(self, ax, title, xlabel, ylabel, xlim, ylim, xmajlocator, xminlocator,
                ymajlocator, yminlocator):
        ax.set_title(f"{title}", fontsize=12)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        if xlabel != None:
            ax.set_xlabel(xlabel, fontsize=10)
        if ylabel != None:
            ax.set_ylabel(ylabel, fontsize=10)    
        if xmajlocator != None: ax.xaxis.set_major_locator(xmajlocator)
        if xminlocator != None: ax.xaxis.set_minor_locator(xminlocator)
        ax.tick_params(axis='both', which='major', labelsize=10)
        if ymajlocator != None: ax.yaxis.set_major_locator(ymajlocator)
        if yminlocator != None: ax.yaxis.set_minor_locator(yminlocator)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

    def style_axes(self, g, xlabel, ylabel, xlim, ylim, xmajmult, xminmult, ymajmult,
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
                self.style_axes_helper(ax, t, xlabel, ylabel, xlim, ylim, xmaj, xmin, ymaj, ymin)
        else:
            self.style_axes_helper(g.figure.axes[0], '', xlabel, ylabel, xlim, ylim, xmaj, xmin, ymaj, ymin)

        g.figure.suptitle(suptitle, fontsize=12, y=0.98, x=0.5)
        g.figure.subplots_adjust(hspace=hspace, wspace=wspace)

    def style_hr(self, g, xlabel, ylabel, xlim, ylim, hspace=0.45, wspace=0.25, 
                suptitle=None, title=None, ymajmult=2, yminmult=1, xmajmult=12):
        self.style_axes(g, xlabel, ylabel, xlim, ylim, xmajmult, 4, ymajmult, yminmult, 
                hspace=hspace, wspace=wspace, suptitle=suptitle, title=title)  

    def style_ms(self, g, xlabel, ylabel, xlim, ylim, xmajmult=500, xminmult=100, 
                ymajmult=2, yminmult=1,hspace=0.45, wspace=0.25, 
                suptitle=None,style=True):
        self.style_axes(g, xlabel, ylabel, xlim, ylim, xmajmult, xminmult, ymajmult, yminmult, hspace=hspace,
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
    
    def style_performance(self, g):
        if g.axes_dict:
            for (name,ax) in g.axes_dict.items(): 
                ax.axhline(y=0, xmin=0, xmax=1, ls="-", color="black", zorder=10, lw=1)
        else:
            g.figure.axes[0].axhline(y=0, xmin=0, xmax=1, ls="-", color="black", zorder=10, lw=1)

    def style_trial(self, g, color='k', alpha=0.3):
        if g.axes_dict: 
            for (name,ax) in g.axes_dict.items():
                ax.fill_between(ax.lines[0].get_data()[0], ax.lines[0].get_data()[1], 
                            color=color, alpha=alpha)
        else:
            ax = g.figure.axes[0]
            ax.fill_between(ax.lines[0].get_data()[0], ax.lines[0].get_data()[1], 
                            color=color, alpha=alpha)
    
    def lineplot(self, data, x=None, y=None, errorbar=None, aspect=1.5, marker='o', 
                mec=None, ms=10, **kwargs):
        
        g = sns.relplot(data, x=x, y=y, kind="line",  aspect=aspect, 
                        errorbar=errorbar, marker=marker, mec=mec,ms=ms,
                        facet_kws={"sharey":False, "sharex":False}, **kwargs)
        return g

    def barplot(self, data, x=None, y=None, hue=None, col=None, color=None, aspect=1,
                palette=None, hue_order=None, legend=None, col_wrap=None, 
                errorbar="se", dodge=False, capsize=0.2, sharex=False, sharey=False,
                native_scale=True, err_kws={"zorder":0.5, "lw":1}, **kwargs):
        
        g = sns.catplot(data, x=x, y=y, palette=palette, col_wrap=col_wrap,
                        hue_order=hue_order, hue=hue, color=color, col=col, 
                        kind="bar", errorbar=errorbar, dodge=dodge, capsize=capsize, 
                        err_kws=err_kws, sharex=sharex, sharey=sharey,aspect=aspect,
                        native_scale=native_scale, legend=legend, **kwargs)
        return g

    def plot_hr(self,data, x="Time (hr)", y="lick", color='k',
                xlim=[-48, 24], ylim=[-10, 10], ymajmult=2, yminmult=1,
                ylabel="", xlabel="Time (hr)", hspace=0.45, 
                wspace=0.25, suptitle='', title=None, **kwargs):
        
        g = self.lineplot(data, x=x, y=y,color=color, **kwargs)
        
        self.style_hr(g, xlabel, ylabel, xlim, ylim, hspace=hspace, wspace=wspace, 
                suptitle=suptitle, title=title,ymajmult=ymajmult, yminmult=yminmult,)  
        
        return g 

    def plot_ant_lickfreq(self, data, x="Time (hr)", y="lick", hue="stimulus", 
                        col="condition", palette=["green", "red"], 
                        hue_order=["stimulus", "blank"],
                        ylim=[0, 12], ylabel="Lick Frequency (Hz)", 
                        **kwargs):
        g = self.plot_hr(data, x=x, y=y, hue=hue, col=col, palette=palette, 
                    hue_order=hue_order, ylim=ylim, ylabel=ylabel, **kwargs)
        
        return g
    
    def plot_ant_perf(self, data, x="Time (hr)", y="lick", col="condition",ylabel=f"Performance\nL_s-L_b", **kwargs):
        g = self.plot_hr(data, x=x, y=y, col=col,ylabel=ylabel, **kwargs) 
        self.style_performance(g)
        return g

    def plot_ant_perf_bar(self, data, x="Time (hr)", y="lick", hue=None, col="condition", 
                        xlim=[-50, 26], ylim=[-10, 10], palette=None, 
                        hue_order=None, color='k', ylabel="Performance", 
                        xlabel="Time (hr)", hspace=0.45, wspace=0.35, legend=True,
                        col_wrap=None, suptitle='', **kwargs):
        
        g = self.barplot(data, x=x, y=y, col=col, col_wrap=col_wrap, hue=hue, 
                    palette=palette, hue_order=hue_order, color=color, 
                    legend=legend, **kwargs)  
        self.style_hr(g, xlabel, ylabel, xlim, ylim, hspace=hspace, wspace=wspace, 
                suptitle=suptitle)
        self.style_performance(g)

        return g

    def plot_lickfreq(self, data, x="Time (ms)", y="lick", hue="stimulus", 
                    col="Time (hr)",color=None,legend=True, col_wrap=6,
                    palette=["green", "red"], hue_order=["stimulus", "blank"], 
                    xlim=[-200, 2000], ylim=[0, 12], ylabel="Lick Frequency",
                    xlabel="Time (ms)", hspace=0.45, wspace=0.25, suptitle=None,
                    xmajmult=500, ymajmult=2, xminmult=100, 
                    yminmult=1, style=True,**kwargs):
        
        g = self.lineplot(data, x=x, y=y, hue=hue, col_wrap=col_wrap, col=col,
                    palette=palette, color=color, hue_order=hue_order, 
                    legend=legend, **kwargs)
        self.style_ms(g, xlabel, ylabel, xlim, ylim, xmajmult=xmajmult, ymajmult=ymajmult,
                xminmult=xminmult, yminmult=yminmult,hspace=hspace, wspace=wspace, 
                suptitle=suptitle,style=style)
        
        return g

    def plot_perf(self, data, x="Time (ms)", y="lick", hue="Time (hr)", 
                col="condition", color='k', legend=True, col_wrap=None,
                palette=None, hue_order=None, xlim=[-200, 2000], ylim=[-10, 10],
                ylabel="Performance", xlabel="Time (ms)", hspace=0.45, 
                wspace=0.25, suptitle=None, ymajmult=2, yminmult=1,xmajmult=500, xminmult=100, **kwargs):
        g = self.lineplot(data, x=x, y=y, hue=hue, col_wrap=col_wrap, col=col,
                    palette=palette, color=color, hue_order=hue_order, 
                    legend=legend, **kwargs)
        self.style_ms(g, xlabel, ylabel, xlim, ylim, hspace=hspace, wspace=wspace, ymajmult=ymajmult,
                suptitle=suptitle, yminmult=yminmult, xmajmult=xmajmult, xminmult=xminmult)
        self.style_performance(g)
        
        return g

    def plot_bar(self,data, x=None, y=None, aspect=1, col=None, 
                    color="#3B3838", palette=None, col_wrap=None, hue_order=None,
                    hue=None, xlim=[-0.5, 2.5], ylim=[0, 600], xlabel=None, hspace=0.45, wspace=0.25,
                    ylabel="", width=0.4, ymajmult=100, yminmult=50,title='', suptitle='',
                    **kwargs):
        g = self.barplot(data, x=x, y=y, aspect=aspect, col=col, color=color, 
                    palette=palette, hue_order=hue_order, hue=hue, 
                    col_wrap=col_wrap, width=width, **kwargs)
        self.style_axes(g, xlabel, ylabel, xlim, ylim, None, None, ymajmult=ymajmult, 
                yminmult=yminmult, hspace=hspace, wspace=wspace, title=title, suptitle=suptitle)

        return g

    def plot_trial_hr(self, data, x="Time (hr)", y="trial no", col="condition", color='k', 
                    hue=None, palette=None, col_wrap=None, hue_order=None, 
                    errorbar='se', err_style='bars', xlabel="Time (hr)", 
                    ylabel="Number of trials", xlim=[-48, 48], ylim=[0, 300], 
                    hspace=0.45, wspace=0.25, suptitle=None, fill=True, title=None, 
                    ymajmult=100, yminmult=50,**kwargs):
        g = self.lineplot(data, x=x, y=y, hue=hue, col=col, color=color, palette=palette, 
                    hue_order=hue_order, errorbar=errorbar, err_style=err_style, 
                    **kwargs)
        self.style_hr(g, xlabel, ylabel, xlim, ylim, ymajmult=ymajmult, yminmult=yminmult,
                hspace=hspace, wspace=wspace, suptitle=suptitle, title=title)
       
        if fill: self.style_trial(g)
        return g

    def connect_lines(self, ax,data, x=None,y=None, hue=None, order=None, hue_order=None):
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

    def plot_bar_strip(self, data, x=None, y=None, row=None, col=None, col_order=None, connect=False,
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
                    self.connect_lines(ax,cdfilt, hue=hue, order=order, x=x, y=y,hue_order=hue_order)
                                
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
                self.connect_lines(ax,data, hue=hue, order=order, x=x,y=y, hue_order=hue_order)
                            
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.tick_params(labelbottom=True, labeltop=False)
            ax.set_xticklabels(ax.get_xticklabels(), **xlab_kws)

        g.tight_layout()
        return g

    def draw_heatmap(self, *args, **kwargs):
        data = kwargs.pop('data')
        yind = kwargs.pop('yind')
        d = data.pivot(index=args[1], columns=args[0], values=args[2])
        for y in yind:
            if y not in d.index:
                d.loc[y] = np.NaN
        d = d.sort_index()
        sns.heatmap(d, **kwargs)

    def row_col_facet(self, data, func_args, func, row, col, row_order=None, col_order=None, **kwargs):
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

    def col_facet(self, data, func_args, func, col, col_wrap=None,col_order=None, **kwargs):
        if col_order:
            cols = col_order
        else:
            cols = data[col].unique()
        if col_wrap:
            if len(cols) % col_wrap == 0:
                nrows = (len(cols) // col_wrap)
            else:
                nrows = (len(cols) // col_wrap) + 1
            fg, axs = plt.subplots(nrows,col_wrap, figsize=(5*col_wrap, 5*nrows))
        else:
            fg, axs = plt.subplots(1, len(cols), figsize=(5*len(cols), 5))
        axs = axs.flatten()
        if col_wrap:
            for i in range(len(cols), nrows*col_wrap):
                fg.delaxes(axs[i])
        fg.tight_layout()
        cax=None
        for j in range(len(cols)):
            plt_data = data[(data[col] == cols[j])]
            cbar = j == (len(cols)-1)
            if cbar:
                cax = axs[len(cols) - 1].inset_axes([1.05, 0, 0.05,1])
            func(*func_args,data=plt_data,ax=axs[j],
                            cbar=cbar,cbar_ax=cax,**kwargs)                            
            axs[j].set_title(f'{cols[j]}')
        fg.tight_layout()
            
        return fg

    def row_facet(self, data, func_args, func, row,row_order=None, **kwargs):
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

    def single_facet(self, data, func_args, func,**kwargs):
        fg, ax = plt.subplots(1, 1, figsize=(5,5))
        cax = ax.inset_axes([1.05, 0, 0.05,1])
        func(*func_args, data=data,ax=ax,
                        cbar=True,cbar_ax=cax,**kwargs) 
        return fg

    def draw_facet(self, data, func_args,  func, row=None, col=None, col_wrap=None, col_order=None, row_order=None,**kwargs):
        yind = data[func_args[1]].unique()
        if row and col:
            cols = col_order if col_order else data[col].unique()
            rows = row_order if row_order else data[row].unique()
            if len(rows) == 1 and len(cols) == 1:
                fg = self.single_facet(data, func_args, func,yind=yind, **kwargs)
                fg.axes[0][0].set_title(f'{rows[0]}, {cols[0]}')
            elif len(rows) == 1:
                fg = self.col_facet(data[data[row] == rows[0]], func_args, func, col, col_wrap=col_wrap, 
                            col_order=col_order,yind=yind,**kwargs)
            elif len(cols) == 1:
                fg = self.row_facet(data[data[col] == cols[0]], func_args, func, row, row_order,yind=yind, **kwargs)
            else:
                fg = self.row_col_facet(data, func_args, func, row, col,
                                row_order=row_order, col_order=col_order, yind=yind,**kwargs)
        elif row:
            rows = data[row].unique()
            if len(rows) == 1:
                fg = self.single_facet(data, func_args, func,yind=yind, **kwargs)
                fg.axes[0].set_title(f'{rows[0]}')
            else:
                fg = self.row_facet(data, func_args, func, row, row_order=row_order,yind=yind,**kwargs)
        elif col:
            cols = data[col].unique()
            if len(cols) == 1:
                fg = self.single_facet(data, func_args, func,yind=yind, **kwargs)
                fg.axes[0].set_title(f'{cols[0]}')
            else:
                fg = self.col_facet(data, func_args, func, col, col_wrap=col_wrap, 
                            col_order=col_order,yind=yind,**kwargs)
        else:
            fg = self.single_facet(data, func_args, func,yind=yind, **kwargs)
        fg.tight_layout()
        return fg

    def plot_pair_strip(self, data, x=None, y=None, row=None, col=None, col_order=None, connect=False,
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
                    self.connect_lines(ax,cdfilt, hue=hue, order=order, x=x, y=y,hue_order=hue_order)                    
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
                    self.connect_lines(ax,cdfilt, hue=hue, order=order, x=x, hue_order=hue_order)                    
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

    def style_retrain(self, g, xlim, ylim, xticklab, xlab, ylab):
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

    def plot_retrain(self, c, means, perf, tot, xlim, xticklab,xlab='Time (hr)', lckylab='Lick Frequency (Hz)', 
                    perfylab="Performance\n$L_s - L_b$", totylab='Number of trials',
                    lckylim=[-0.1, 10], perfylim=[-6,6], totylim=[0, 400], aspect=2, **kwargs):
        g1 = None
        if not means.empty:
            g1 = self.plot_ant_lickfreq(means[(means['cond'] == c) & (means['test_time'] == 'train')],
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
            self.style_retrain(g1, xlim, lckylim, xticklab, xlab, lckylab)
        g2 = None
        if not perf.empty:
            g2 = self.plot_ant_perf(perf[(perf['cond'] == c) & (perf['test_time'] == 'train')],
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
            self.style_retrain(g2, xlim, perfylim, xticklab,xlab, perfylab)
        g3 = None
        if not tot.empty:
            g3 = self.plot_trial_hr(tot[(tot['cond'] == c) & (tot['test_time'] == 'train')],
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
            self.style_retrain(g3, xlim, totylim, xticklab, xlab, totylab)
            for i in range(len(g3.axes[0][0].lines)):
                if i % 2 == 0:
                        g3.axes[0][0].fill_between(g3.axes[0][0].lines[i].get_data()[0], g3.axes[0][0].lines[i].get_data()[1], 
                                                color='k', alpha=0.3)
        return (g1, g2, g3)
    
    def generate_ant_lf(self):
        if self.antlf:
            try:
                ant_means = self.data['fixed means filtered']
            except KeyError:
                errdisp = tk.Toplevel()
                errlab = tk.Label(errdisp, text='No fixed analysis')
                errlab.grid()
                return
            ant_means['Time (hr)'] = ant_means['Time (hr)'] + 2
            xlim = [ant_means['Time (hr)'].min() - 2, ant_means['Time (hr)'].max() + 2]
            g = self.plot_ant_lickfreq(ant_means, xlim=xlim, aspect=2, ms=7)
            g.figure.subplots_adjust(left=0.05, bottom=0.05, top=0.95, right=0.95)
            if self.disp:
                self.lfcan = MplEmbedPlot(g.figure, master=self.plt_container)
                self.plt_container.add(self.lfcan.frame, text='Fixed Lick Frequency')
            if self.save:
                g.savefig(self.savedir+"\\fixed_lickfreq.png", transparent=True)

    def generate_ant_perf(self):
        if self.antperf:
            try:
                perf = self.data['fixed performance filtered']
            except KeyError:
                errdisp = tk.Toplevel()
                errlab = tk.Label(errdisp, text='Fixed Performance: No fixed analysis')
                errlab.grid()
                return
            perf['Time (hr)'] = perf['Time (hr)'] + 2
            xlim = [perf['Time (hr)'].min() - 2, perf['Time (hr)'].max() + 2]

            g = self.plot_ant_perf(perf, ylabel="Performance", aspect=2, xlim=xlim,ms=7, ylim=[-10, 10])
            g.figure.subplots_adjust(left=0.05, bottom=0.05, top=0.95, right=0.95)
            if self.disp:
                self.perfcan = MplEmbedPlot(g.figure, master=self.plt_container)
                self.plt_container.add(self.perfcan.frame, text='Fixed Performance')
            if self.save:
                g.savefig(self.savedir+"\\fixed_perf.png", transparent=True)

    def genereate_fullckhr(self):
        if self.instlf:
            try:
                means_full = self.data['inst means filtered']
            except KeyError:
                errdisp = tk.Toplevel()
                errlab = tk.Label(errdisp, text='No instantaneous analysis')
                errlab.grid()
                return
            self.fllf_frame = ttk.Frame(self.plt_container)
            self.plt_container.add(self.fllf_frame, text='Intstantaneous Lick Freq')
            self.fllf_nte = ttk.Notebook(self.fllf_frame)
            self.fllfplots = dict()
            self.fllf_nte.pack()
            for c in means_full['condition'].unique():
                g = self.plot_lickfreq(means_full[means_full['condition'] == c], suptitle=c,
                                xlim=[-700, 2000], col="Time (hr)",col_wrap=6, ms=1, lw=0.7, 
                                ymajmult=5, yminmult=2.5, xminmult=250, wspace=0.2, hspace=0.2)
                g.tick_params(labelsize=8, pad=0.5)
                i=0
                for n, ax in g.axes_dict.items():
                    ax.set_xlabel(ax.get_xlabel(), fontsize=8)
                    ax.set_ylabel(ax.get_ylabel(), fontsize=8)
                    ax.set_title(f'Time (hr):{ax.get_title()}', fontsize=8)
                    if i % 6 != 0:
                        ax.set_ylabel('')
                    if i < len(g.axes)-6:
                        ax.set_xlabel('')
                    ax.tick_params(axis='x', which='major', labelrotation=45)
                    i +=1
                nrow = math.ceil(len(g.axes)/6)
                g.figure.set_size_inches(6*1.5+1, nrow*1.5)
                g.legend.set_in_layout(True)
                g.figure.tight_layout()

                g.legend.set_title("Stimulus Type", prop=mpl.font_manager.FontProperties(size=8))
                for t in g.legend.get_texts():
                    t.set_size(8)
                g.legend.set(loc='center right',bbox_to_anchor=(1.01, 0.5))
                g.figure.subplots_adjust(left=0.05, bottom=0.15,right=0.9)
                
                if self.disp:
                    lfdisp = MplEmbedPlot(g.figure, master=self.fllf_nte)
                    self.fllfplots[c] = lfdisp
                    self.fllf_nte.add(lfdisp.frame, text=f'Instentaneous Lick Freq {c}')
                if self.save:
                    g.savefig(self.savedir.get()+f"\\lckfreqfull_{c}.png", transparent=True)

    def generate_fullpfhr(self):
        if self.instperf:
            try:
                perf = self.data['inst performance filtered']
            except KeyError:
                errdisp = tk.Toplevel()
                errlab = tk.Label(errdisp, text='No instantaneous analysis')
                errlab.grid()
                return
            flperf = ttk.Frame(self.plt_container)
            self.plt_container.add(flperf, text='Intstantaneous Performance')
            self.flperf_nte = ttk.Notebook(flperf)
            self.flperfplots = dict()
            self.flperf_nte.pack()
            for c in perf['condition'].unique():
                g = self.plot_perf(perf[perf['condition'] == c], xlim=[-700, 2000], 
                            col="Time (hr)", col_wrap=6, hue=None, color='k',suptitle=c,
                            ms=1, lw=0.7, ymajmult=5, yminmult=2.5, xminmult=250, wspace=0.2,hspace=0.2 )
                g.tick_params(labelsize=8, pad=0.5)
                i=0
                nrow = math.ceil(len(g.axes)/6)
                for n, ax in g.axes_dict.items():
                    ax.set_xlabel(ax.get_xlabel(), fontsize=8)
                    ax.set_ylabel(ax.get_ylabel(), fontsize=8)
                    ax.set_title(f'Time (hr):{ax.get_title()}', fontsize=8)
                    if i % 6 != 0:
                        ax.set_ylabel('')
                    if i < len(g.axes)-6:
                        ax.set_xlabel('')
                    ax.tick_params(which='major', axis='x', labelrotation=45)
                    i +=1
                g.figure.set_size_inches(6*1.5, nrow*1.5)
                g.tight_layout()
                g.figure.subplots_adjust(left=0.05, bottom=0.15)
                if self.disp:
                    flperfcan = MplEmbedPlot(g.figure, master=self.flperf_nte)
                    self.flperfplots[c] = flperfcan
                    self.flperf_nte.add(flperfcan.frame, text=f'Instentaneous Performance {c}')
                if self.save:
                    g.savefig(self.savedir +f"\\perffull_{c}.png", transparent=True)

    def generate_htmplf(self):
        if self.htmplf:
            try:
                self.data['inst means filtered']
            except KeyError:
                errdisp = tk.Toplevel()
                errlab = tk.Label(errdisp, text='No instantaneous analysis')
                errlab.grid()
                return
            gpcols = ['condition', 'Time (hr)', 'Time (ms)','stimulus']
            means_full_means = self.get_means(self.data['inst means filtered'], gpcols)
            htmp_data = means_full_means.sort_values(['Time (hr)', 'Time (ms)'])
            htmp_data = htmp_data[(htmp_data['Time (ms)'] < 2000) & (htmp_data['Time (ms)'] > -700)]
            htmp_args=['Time (ms)', 'Time (hr)', 'lick']
            fg = self.draw_facet(htmp_data, htmp_args, self.draw_heatmap, col='stimulus', row='condition',
                    col_order=['stimulus', 'blank'],square=True,
                    cmap='viridis', vmin=0, vmax=10, xticklabels=3)
            fg.suptitle('Instantaneous licking', y=0.98)
            fg.subplots_adjust(left=0.05, bottom=0.15, top=0.9, right=0.95, wspace=0, hspace=0.75)
            if self.disp:
                self.fdisp = MplEmbedPlot(fg.figure, master=self.plt_container)
                self.plt_container.add(self.fdisp.frame, text='Inst Lick Freq Heatmap')
            d = htmp_data.pivot(index=["condition",  'stimulus',"Time (hr)"], columns="Time (ms)", values='lick')
            if self.save:
                fg.savefig(self.savedir+"\\lckfreqheatmap.png", transparent=True)
                d.to_csv(self.savedir+"\\lckfreqheatmap.csv")
            else:
                if not os.path.isdir(self.loaddir + "\\heatmap"):
                    os.mkdir(self.loaddir + "\\heatmap")
                d.to_csv(self.loaddir+"\\heatmap\\lckfreqheatmap.csv")

    def generate_htmppf(self):
        if self.htmppf:
            try:
                self.data['inst performance filtered']
            except KeyError:
                errdisp = tk.Toplevel()
                errlab = tk.Label(errdisp, text='No instantaneous analysis')
                errlab.grid()
                return
            gpcols = ['condition', 'Time (hr)', 'Time (ms)']
            perf_full_means = self.get_means(self.data['inst performance filtered'], gpcols)
            htmp_data = perf_full_means.sort_values(['Time (hr)', 'Time (ms)'])
            htmp_data = htmp_data[(htmp_data['Time (ms)'] < 2000) & (htmp_data['Time (ms)'] > -700)]
            htmp_args=['Time (ms)', 'Time (hr)', 'lick']
            fg = self.draw_facet(htmp_data, htmp_args, self.draw_heatmap, col='condition',square=True,
                    cmap='coolwarm', vmin=-6, vmax=6, col_wrap=6, xticklabels=3)
            fg.suptitle('Instantaneous performance', y=0.98)
            fg.subplots_adjust(left=0.05, bottom=0.1, top=0.9, right=0.95, wspace=0.25, hspace=0.25)
            if self.disp:
                self.htmppfdisp = MplEmbedPlot(fg.figure, master=self.plt_container)
                self.plt_container.add(self.htmppfdisp.frame, text='Inst Perf Heatmap')
            if self.save:
                fg.savefig(self.savedir+"\\perfheatmap.png", transparent=True)
 
    def generate_nthpart(self):
        if self.nthpart:
            try:
                last20 = self.data['part means']
            except KeyError:
                errdisp = tk.Toplevel()
                errlab = tk.Label(errdisp, text='No part analysis')
                errlab.grid()
                return
            last20["Day"] = last20["Day"].apply(self.day_to_label)
            g = self.plot_pair_strip(last20, y="lick",x="stimulus", order=["stimulus", "blank"],
                            hue="stimulus", palette=["green", "red"], ylabel='Lick frequency (Hz)',
                            hue_order=["stimulus", "blank"],
                            row="condition", col="Day",jitter=False, s=7, aspect=0.5, legend=False, connect=True)
            g.figure.subplots_adjust(left=0.05, bottom=0.05, top=0.95, right=0.95)
            for (r, c), ax in g.axes_dict.items():
                ax.set_title(f'{r} {c}', fontsize=10)
            if self.disp:
                self.nthdisp = MplEmbedPlot(g.figure, master=self.plt_container)
                self.plt_container.add(self.nthdisp.frame, text='Nth Part')
            if self.save:
                g.savefig(self.savedir+"\\nthpart.png", transparent=True)

    def generate_trialshr(self):
        if self.trialshour:
            try:
                trials = self.data['fixed trial counts']
            except KeyError:
                errdisp = tk.Toplevel()
                errlab = tk.Label(errdisp, text='No fixed analysis')
                errlab.grid()
                return
            tots = self.sum_trials(trials, ['condition', 'animal', 'Time (hr)'],[], "trial no")
            xlim = [tots['Time (hr)'].min() - 2, tots['Time (hr)'].max() + 2]
            ymax = math.ceil(tots['trial no'].max()/100.0)*100
            g = self.plot_trial_hr(tots, xlim=xlim,ylim=[0, ymax], title=None, ms=7, lw=1,err_kws={"lw":0.5}, ylabel='Trial no')
            g.figure.subplots_adjust(left=0.05, bottom=0.05, top=0.95, right=0.95)
            if self.disp:
                self.trialsdisp = MplEmbedPlot(g.figure, master=self.plt_container)
                self.plt_container.add(self.trialsdisp.frame, text='Trials by hour')
            if self.save:
                g.savefig(self.savedir +"\\trialsbyhour.png", transparent=True)
    
    def day_to_label(self,day):
        day = int(day)
        if day < 0:
            return "ACC " + str(day + 1)
        return "SAT " + str(day + 1)
    
    def sum_trials(self,data, group, keep, key):
        '''Get total number of trials grouped by given group. Only use on aggregated
        trial data.
        
        Parameters:
        data --- data frame containing trial data. Must contain all columns in
                `group`, `keep`, and `key`.
        group --- list of columns to group by
        keep --- list of columns to keep
        key --- column containing trial counts'''
        tot = data.groupby(group)[key].sum().rename('trial no')
        keep = data.groupby(group).first()[keep]
        res = pd.concat([tot, keep], axis=1).reset_index()
        return res
    
### Main Analysis Code ###
class RunAnalysisWindow(tk.Toplevel):
    def __init__(self, master=None, **kwargs):
        self.tl = tk.Toplevel.__init__(self,**kwargs)
        self.title('Analysis Progress')
        self.geometry('500x500+100+100')
        
        self.mssg = StringVar(value="Analysis\n--------")
        self.proglab = ttk.Label(self, textvariable=self.mssg)
        self.proglab.grid(row=0, column=0, padx=2, pady=2)

        self.in_fold = master.files_var.get()
        self.metadata_n = master.metadata_var.get()
        print(self.metadata_n)
        self.metadata = pd.read_excel(self.metadata_n)
        self.time_bin = int(master.time_bin.get()) 
        self.keep = master.keep_var.get().split(', ')
        
        if len(master.values_var.get().split(', ')) == 1:
            self.values = [master.values_var.get()]
        else:
            self.values = master.values_var.get().split(', ') 

        self.ifinst = master.inst_anal.get() 
        self.iffixed = master.fixed_anal.get()
        self.fixedstart = int(master.start_val.get())
        self.fixedend = int(master.end_val.get())
        self.ifnthfrac = master.nth_frac_val.get()
        self.num, self.denom = int(master.num_val.get().split('/')[0]), int(master.num_val.get().split('/')[1])
        self.min_trials = int(master.min_var.get())
        self.min_blank = int(master.min_blank_var.get()) 
        self.min_stim = int(master.min_stimulus_var.get())
        self.add = master.add_val.get()
        self.out_fold = master.add_file_label_var.get()
        self.fold_type = master.type_var.get() 
        self.as_cond = (self.fold_type=="Condition")

        self.mins = {
            "min_trials": self.min_trials,
            "min_blank": self.min_blank,
            "min_stimulus": self.min_stim
            }
        self.index = ["condition", "animal"]
        self.index_trial = ['condition', 'animal', 'trial no']
        self.full_keys = {
                "millis bin":"puff delta",
                "Time (ms)":"Time (ms)",
                "stimulus type":"stimulus",
                "water delivery":"water",
                "trial no":"trial no",
                "time bin":"delta",
                "licks":"lick",
                "timestamp":"timestamp",
                "delay":"delay",
                "acc":"acc",
                'offset':'offset'
            }
        
        self.full_suffix = ['data', 'means', 'means_filtered', 
                    'performance', 'performance_filtered', 
                    ]
        self.ant_suffix = self.full_suffix + ['trial counts', 'trials summed day'] 

        self.last20_suffix = ['data', 'means', 'performance', 'trial counts']
        
        self.out_name = master.name_var.get()
        self.cond_name = Path(self.in_fold).parts[-1]
        if self.out_name == '':
            self.out_name=self.cond_name    

        if self.add:
            analysis_output = os.listdir(self.out_fold)
            analysis_output.sort()
            self.full_out = [self.out_fold + '\\' +  analysis_output[i] for i in range(7,12)]
            self.ant_out = [self.out_fold + '\\' +  analysis_output[i] for i in range(0,7)]
            self.last20_out = [self.out_fold + '\\' +  analysis_output[i] for i in range(12,16)]
        else:
            self.full_out = [self.out_fold + '\\inst_'  + self.out_name + '_' + suff + '.csv' for suff in self.full_suffix]
            self.ant_out = [self.out_fold + f'\\fixed {self.fixedstart}-{self.fixedend}_'  + self.out_name + '_' + suff + '.csv' for suff in self.ant_suffix]
            self.last20_out = [self.out_fold + f'\\{self.num}-{self.denom} part _'  + self.out_name + '_' + suff  + '.csv' for suff in self.last20_suffix]
        self.out_files = self.full_out + self.ant_out + self.last20_out

        self.plot = master.analyplot.get()


        
        self.analysis()
    
    def analysis(self):
        self.after(1, lambda : self.run_formatter())    

    #######################################
    #  Functions to Load and Format Data  #
    #######################################
    def time_from_file(self, f):
        "Convert filename to timestamp."
        dt = re.findall("\d\d_\d\d_\d\d_?~?T_?~?\d\d_\d\d_\d\d", f)[0]
        dt = dt.replace("_", "-", 2).replace("~", "-", 2).replace("_T_", " ").replace("-T-", " ").replace("_", ":")
        return datetime.datetime.strptime(dt,'%m-%d-%y %H:%M:%S')

    def set_noon(self,t):
        "Get noon on start of training as timestamp"
        return datetime.datetime(t.year, t.month, t.day, 12, 0, 0)

    def load_file(self, file_path):
        '''Loads, formats, and returns lick frequecny data as a data frame.
        
        Loads data from csv format. Extracts the start time from the filename and uses it to format timestamps from milliseconds since start of file to datetimes.
        Converts delays at beginning of trial to timedeltas.

        Parameters:
        filename --- path of file to load (string) '''

        data = pd.read_csv(file_path, header=None)

        # rename columns
        mp = {0:"timestamp",1:"poke", 2:"lick", 3:"condition code", 4:"delay"}
        if len(data.columns) == 6:
            mp = {0:"timestamp",1:"poke", 2:"lick", 3:"condition code", 4:"delay", 5:"stimulus"}
        data = data.rename(columns=mp)
        
        # convert time columns to correct type
        dt = self.time_from_file(file_path)
        # convert to int in ms to avoid floating point rounding errors in datetime conversion
        data['timestamp'] = (data['timestamp']*1000).astype(int)
        data["timestamp"] = pd.to_datetime(data["timestamp"], unit="ms", origin=dt)
        data["delay"] = pd.to_timedelta(data["delay"], unit="ms")

        # get noon for timebin alignment
        data['offset'] = self.set_noon(dt)

        return data

    def get_trials(self,trial):
        '''Generates a pandas Series of x of length len.
        
        Used to label all samples in a trial with the trial number.

        Parameters
        trial --- (value, length)'''
        (x, len) = trial
        # generates a list of x of length len
        return pd.Series([x for i in range(len)], dtype="int")

    def enumerate_trials(self,data, start=0):
        '''Labels each sample in a trial with the trial type and the trial number, numbering trials consecutively from `start`.
        The number of trials, and the number of each trial type, matches v16 of the matlab analysis. Returns the input dataframe with labeled trials.
        '''
        # get first sample in each trial
        # first sample, all samples where current water code is not 7 (timeout) and previous is 7,  last sample
        trial_boundary = pd.concat([data.iloc[[0]], data[(data["condition code"]!=7) & (data["condition code"].shift() == 7)], data.iloc[[data.shape[0]-1]]])
        # indexes of first sample in each trial
        trial_boundary_indicies = pd.Series(trial_boundary.index)

        # get number of samples per trial 
        # (difference in index number between previous sample and current sample)
        num_samples = trial_boundary_indicies.diff().fillna(0).astype('int').tolist()
        # label with number of samples with index
        trial_count = pd.Series(enumerate(num_samples))

        # label all samples in a trial with its trial number
        trial_count = trial_count.map(self.get_trials).explode()[1:]
        trial_count.index = range(0, trial_count.shape[0])
        trial_count = trial_count + start
        
        # add trial labels to data
        data["trial no"] = trial_count
        # set trial no for last sample
        data.loc[data.shape[0] - 1, "trial no",] = trial_boundary.shape[0] - 1 + start

        return data

    def label_trials_stimulus(self,data):
        '''Label samples as stimulus or blank.'''
        data["stimulus"] = data["stimulus"].replace({0:"blank", 1:"stimulus"})
        return data

    def label_trials_water(self,data):
        '''Return data with all trials labeled as water or no water.
        
        Trial types are decided by the code in the 3rd column (if the trial contains a 3, it is a water trial). 
        Requires trials to already be grouped individually and labeled with trial number.
        '''
        # get all trials labeled with a 3 (water trials)
        go = data.groupby(["trial no"]).filter(lambda x: (x["condition code"]==3).any())
        # label all these trials as water
        go["water"] = ["water" for i in range(go.shape[0])]
        data["water"] = go["water"]
        # label remaning trials as blank
        data["water"] = data["water"].fillna(value="no water")
        return data

    def label_trials(self,data):
        '''Label trials with both stimulus or blank and water or no water distinction.'''
        data = self.label_trials_water(data)
        if "stimulus" in data.columns:
            data = label_trials_stimulus(data)
        else:
            data["stimulus"] = data["water"]
            data["stimulus"] = data["stimulus"].replace({"water":"stimulus", "no water":"blank"})
        data["type"] = data["water"] + " & " + data["stimulus"]
        return data

    def drop_last_sample(self,data):
        '''Drops last sample of each trial. This has the same timestamp as 
        the second-to-last sample and only indicates the trial is over.'''
        data = data[data.groupby(["animal", "trial no"]).cumcount(ascending=False) > 0]
        return data

    def make_animal_df(self,andir, animal_name):
        '''Load and format all data files for one animal and return as a data frame.

        Requires:
        --- `andir` only contains data files (may contain directories but not other files)
        --- metadata file is provided
        --- metadata contains animal ID in column named "Animal ID"
        --- metadata contains acclimation time in column named "acc"

        If name of animal in metadata file does not match `animal_name` no data will
        be loaded. Column names in addition to "acc" can be provided to load metadata 
        from those columns in the order they are provided. Prints a warning message
        for all additional columns that were provided but had no data, and returns
        those columns as NaN. For animals that have training split into multiple files,
        the entire filename except for the timestamp _must_ be identical, or they
        will not be loaded in the correct order.

        Parameters:
        andir --- directory containing the data files for that animal (and no other files)
        metadata --- dataframe containing metadata about an animal associated with animal id
                    if animal id column is not named exactly 'Animal ID', will fail
        animal_name --- ID code of animal
        meta_vals --- list of columns to include metadata from. Must include 'acc' for length of acclimation (in days).'''
        
        if self.metadata.empty:
            print("Must include metadata file including at least acclimation time")
            mssg.set("Must include metadata file including at least acclimation time")
            return
        if 'acc' not in self.metadata.columns:
            print("Column indicating acclimation time must be named 'acc'")
            self.mssg.set(self.mssg.get() + "\nColumn indicating acclimation time must be named 'acc'")
            return
        if 'acc' not in self.keep:
            print("Must include acc in list of columns to be selected from metadata")
            self.mssg.set(self.mssg.get() + "\nMust include acc in list of columns to be selected from metadata")
            return

        fs = os.listdir(andir)
        # ensure files concatenated in time order
        # only works if rest of filename is identical
        fs.sort()
        animal = []        
        # load all files in an animal's folder
        start = 0
        for f in fs:
            f_path = andir + "\\" + f
            if (os.path.isfile(f_path)): 
                try:     
                    data = self.load_file(f_path)
                except UnicodeDecodeError:
                    print(f"Check that only behavior files are in {andir}")
                    self.mssg.set(self.mssg.get() + "\nCheck that only behavior files are in {andir}")
                    return
                # labeling trials
                data = self.enumerate_trials(data, start)
                start = data.tail(1)["trial no"].reset_index(drop=True)[0]
                data = self.label_trials(data)
                animal.append(data)
        # load metadata
        if animal == []:
            print(f"{animal_name} has no behavior files")
            self.mssg.set(self.mssg.get() + f"\n{animal_name} has no behavior files. Check if this folder should be run as an animal or as a condition.")
            return pd.DataFrame()
        else:
            animal = pd.concat(animal, ignore_index=True)
            try:
                self.metadata[self.metadata["Animal ID"] == animal_name].empty
            except KeyError:
                print("Animal ID column must be named 'Animal ID'")
                self.mssg.set(self.mssg.get() + "\nAnimal ID column must be named 'Animal ID'")
                return
            if not (self.metadata[self.metadata["Animal ID"] == animal_name].empty):
                an_meta =  self.metadata[self.metadata["Animal ID"] == animal_name].reset_index()
                fail_keys = []
                for key in self.keep:
                    try:
                        animal[key] = an_meta.loc[0, key]
                        if an_meta.loc[:,key].isna().any():
                            if key == 'acc':
                                print(f'Acclimation time (required) is missing for {animal_name}. Skipping animal.')
                                self.mssg.set(self.mssg.get() + f'\nAcclimation time (required) is missing for {animal_name}. Skipping animal.')
                                return pd.DataFrame()
                            else:
                                fail_keys += [key]
                    except KeyError:
                        fail_keys += [key]
                        if key == 'acc':
                            print("Metadata file must include acclimation time in column named 'acc'.")
                            self.mssg.set(self.mssg.get() + "\nMetadata file must include acclimation time in column named 'acc'.")
                            return
                if len(fail_keys) > 0:
                    print("%s %s not in metadata file" % (animal_name, ', '.join(fail_keys)))     
                    self.mssg.set(self.mssg.get() + ("\n%s %s not in metadata file" % (animal_name, ', '.join(fail_keys)))) 
            else:
                print(f"{animal_name}: no metadata - must have at least acclimation time. Skipping animal.")
                self.mssg.set(self.mssg.get() + f"\n{animal_name}: no metadata - must have at least acclimation time. Skipping animal.")
                return pd.DataFrame()
            
            animal["animal"] = animal_name
            animal['offset'] = animal.loc[0, 'offset']
            animal[fail_keys] = pd.NA
            # one-hot encode lick frequency
            animal.loc[animal["lick"] == 2,"lick"] = 1
            # convert acclimation time to timedelta
            tmp = pd.to_timedelta(animal["acc"], unit="day")
            animal.drop(columns="acc")
            animal["acc"] = tmp
            animal = self.drop_last_sample(animal)
            return animal
            
    #######################################
    #      Analysis Helper Functions      #
    #######################################

    def time_to_float(self,data, name, key, frequency):
        '''Return dataframe with a new column representing `key` column as a float.
        
        Parameters:
        data --- dataframe containing trial data. Must contain `key`.
        name --- name of new column
        key --- name of timedelta64 column to convert
        frequency --- string frequency at which to convert timedelta (eg. ms, s, h)'''
        data[name] = data[key].to_numpy(dtype=f"timedelta64[{frequency}]").astype("float")
        return data

    def puff_delta(self,data, index, key_time, key_delay, name):
        '''Returns dataframe with new column of the given name with caluculated 
        time from each sample in a trial to airpuff delivery in that trial.
        
        Assumes SAT trial structure: IR beam break, then random delay between 200 
        and 800 ms, then 500ms of air puff delivery, then 500ms fixed delay, then 
        water delivery and timeout for 2000 ms or until animal stops licking 
        before next trial. All samples during the random delay period have a non-zero delay. 
        The delay value is taken from the first sample in a trial, which is 
        guaranteed to be during this random delay period. Start of air delivery
        is calculated as the time of the first sample in a trial plus the length of 
        the random delay. This value is subtracted from all samples in the trial to 
        get time from each sample to airpuff delivery.

        Parameters:
        data --- data frame containing trial data. Must contain `key_time`, `key_delay`,
                and all columns in `index`.
        index --- list of columns to group data by
        key_time --- name of column containg sample time
        key_delay --- name of column containing delay time for that trial
        name --- name of column to write delta to
        '''

        data = data.set_index(index)
        grouped = data.groupby(index)
        data["air_start"] = grouped[key_time].first() + grouped[key_delay].first()
        data[name] = data[key_time] - data["air_start"]
        data = data.drop(columns="air_start")

        return data.reset_index()

    def get_first_sample(self,data, index, key, name):
        '''Returns dataframe with new column that has first row of grouped data 
        iterated across all rows in the group. 

        Used to copy timestamp at trial initation across all samples in a trial.

        Parameters:
        data --- data frame containing trial data. Must contain `key` and all 
                columns in `index`.
        index --- list of columns to group data by
        key --- column to take start time value from
        name --- name of column to write first sample to'''
        
        data = data.set_index(index)
        data[name] = data.groupby(index)[key].first()
        data = data.reset_index()
        return data

    def fixed_window_get_lickfreq(self,data,  index, values_mean, values_first, puff, start, end):
        '''Filter data to a fixed window of trial and return lick frequency in that
        window. Lick frequency (in Hz) is calculated as the total number of licks in
        the fixed window divided by the length of the window.

        Parameters:
        data --- data frame containing trial data. Must contain all columns in
                `values`, `keep`, and `index`
        values_mean --- list of columns to calculate lick frequeny
        values_first --- list of columns to keep
        index --- list of columns to index and group by
        puff --- name of column containing time to airpuff delivery
        freq_window --- length of window to calculate rolling window with in ms
        start --- start of fixed window relative to airpuff delivery in ms
        end --- end of fixed window relative to airpuff delivery in ms
        '''
        freq_window = end - start
        data = data[data[puff] >= pd.to_timedelta(start, unit="ms")]
        data = data[data[puff] <= pd.to_timedelta(end, unit="ms")]
        group = data.groupby(index)
        licks = group[values_mean].sum()/(freq_window/1000)
        keep = group[values_first].first()
        data_freq = pd.concat([licks, keep], axis=1).reset_index()
        return data_freq

    def bin_by_time(self,data, bin_size, bin_unit, index, values_agg, values_first, key, 
                    origin="start_day",offset=None, fn="mean"):
        '''Groups trials by time in given size and unit and apply provided function, and returns data with key 
        column replaced with time bin.

        If `fn` is "lickfreq", calculates lickfreq (number of licks/timebin size) for
        each timebin. If fn is `mean`, `std`, `var`, or `se`, calculates the respective
        measure. All other calculations must be a user-provided function.

        Parameters:
        data --- data frame containing trial data. Must contain all columns in 
                `index`, `values_mean`, `values_first` and `key`.
        bin_size --- size of bin to group by, in bin_unit (int)
        bin_unit --- unit of bin to group by (e.g. ms, min, hr)
        index --- list of columns to group data by
        values_agg --- set of columns to aggregate values for
        values_first --- set of columns to keep
        key --- column to group with, replaced with time bin
        origin --- start point to calculate timebins from
        offset --- amount to offset timebins from origin'''

        gp = pd.Grouper(key=key, freq=pd.to_timedelta(bin_size, unit=bin_unit), 
                        offset=offset, origin=origin)
        # append resampler to index group
        index = index.copy()
        index.append(gp)
        group = data.groupby(index)
        # calculate lick frequency
        if fn == 'lickfreq':
            licks = group[values_agg].sum()/(bin_size/1000)
        # handle some other common aggregations that are faster with builtins
        elif fn == 'mean':
            licks = group[values_agg].mean()
        elif fn == 'std':
            licks = group[values_agg].std()
        elif fn == 'var':
            licks = group[values_agg].var()
        elif fn == 'se':
            licks = group[values_agg].se()
        # all other calcualtions
        else:
            licks = group[values_agg].apply(fn)    

        keep = group[values_first].first()
        data_sample = pd.concat([licks, keep], axis=1).reset_index()
        return data_sample

    def delta(self,data, index, key, key_acc, key_offset, name):
        '''Calculates time from start of SAT to trial and returns data in given 
        column.

        For the given group, subtracts the values of key_time and key_acc
        in the first row of the group from all values of key_time in the group.
        If used to calculate time from start of SAT, must provide the timebin each
        row is in as the key.

        Parameters:
        data --- data frame containing trial data
        index --- list of columns to index data
        gp --- column or list of columns to group by
        key --- column with bin times to use for calculation
        name --- name of column to write delta to'''

        data = data.set_index(index)
        delta = data[key] - (data[key_offset] + data[key_acc])
        data[name] = delta
        return data.reset_index()

    def mean_bin(self,data, index, value, keep):
        '''Return mean of columns in `value` grouped by `index`, with columns in 
        `keep` retained.
        
        Parameters:
        data --- data frame containing trial data. Must contain columns in `index`,
                `value`, and `keep`.
        index --- list of columns to group by
        value --- list of columns to average
        keep --- list of columns to retain
        '''

        data = data.groupby(index)
        lick = data[value].mean()
        keep = data[keep].first()

        data = pd.concat([lick, keep], axis=1).reset_index()
        return data

    def performance(self,data, index, keep, key_trial, key_value, cond0, cond1):
        '''Calculate performance (difference between trial types)
        and return as new data frame.

        Groups data by index and given trial type column, then computes 
        cond1 - cond0 for each group to get performance. Takes the first row from 
        cond1 group for kept columns.
        
        Parameters:
        data --- data frame containing trial data. Must contain all columns in 
                `index`, `keep`, `key_trial`, and `key_value`.
        index --- list of columns to group data by
        keep --- list of columns to keep
        key_trial --- column containing trial type labels
        key_value --- column containing data to calculate performance
        cond0 --- name of first trial type
        cond1 --- name of second trial type'''

        group = data.set_index(index).groupby(key_trial)

        blank = group.get_group(cond0)[key_value]
        stim = group.get_group(cond1)[key_value]

        data_perf = stim - blank
        keep = group.get_group(cond1)[keep]
        data_perf = pd.concat([data_perf, keep], axis=1)

        return data_perf.reset_index()

    def trial_counts(self,data, index, keep, key):
        '''Returns total number of rows in key for data grouped by index. 
        
        Parameters:
        data --- data frame containing trial data. Must contain all columns in 
                `index`, `keep`, and `key`.
        index --- list of columns to group data by
        keep --- list of columns to keep
        key --- column to calculate trial counts from
        '''

        gp = data.groupby(index)
        counts = gp[key].nunique()
        keep = gp.first()[keep]

        data = pd.concat([counts, keep], axis=1).reset_index()
        return data

    def get_nth_part_trials(self,start, end, gp):
        '''Get set of trials from gp starting at `start` and ending at `end`. For 
        getting set of trials which represent nth fraction of a day.

        Parameters:
        gp --- group to get trials for
        start --- dataframe of starting trial no correpsonding to group names
        end --- datafram of ending trial no corresponding to group names
        '''
        if ~gp.empty:
            return gp[(gp["trial no"] > start.loc[gp.name]) & (gp["trial no"] <= end.loc[gp.name])]
        else:
            return []

    def get_nth_part_day(self,data, frac, nth_part):
        '''Get set of trials representing the nth fraction of all trials an animal
        did in a day.
        
        Parameters:
        data --- data frame containing trial data
        frac --- denomenator of fraction of trials
        nth_part --- numenator of fraction of trials'''
    # TODO: doesn't work correctly with string days -- FIXED?
        if nth_part > frac:
            print("Partition to return must be less than the number of partitions")
            return
        
        # total number of trials per day
        total_trials_day= data.sort_values(['condition', 'animal', 'trial no']).groupby(['condition', 'animal', 'Day']).last()['trial no']
        total_trials_day1 = total_trials_day.reset_index().groupby(['condition', 'animal']).first().reset_index().set_index(['condition', 'animal', 'Day'])['trial no']
        total_trials_day = total_trials_day.groupby(['condition', 'animal']).diff().fillna(total_trials_day1)
        
        # number of trials in a partition per day, rounded down
        trials_per_part = np.floor(total_trials_day*(1/frac))
        
        # first and last trial in partition for each day
        start = (trials_per_part*(nth_part-1))
        end = (trials_per_part*nth_part)
        
        # adjust first and last trial numbers for number of trials done on previous 
        # days
        for i in range(1, len(data["Day"].unique()) + 1):
            start = start + total_trials_day.groupby("animal").shift(i).fillna(0)
            end = end + total_trials_day.groupby("animal").shift(i).fillna(0)

        # convert to int to allow indexing
        start = start.astype(int)
        end = end.astype(int)
        
        # get trials in partition
        part = data.set_index(['condition', "animal", "Day"]).groupby(['condition', "animal", "Day"], group_keys=False).apply(lambda x: self.get_nth_part_trials(start, end, x), include_groups=False)
        return part.reset_index()

    def day_to_label(self,day):
        '''Convert day as a float time to start of SAT representation to a string representation.'''
        day = int(day)
        if day < 0:
            return "ACC " + str(day + 1)
        return "SAT " + str(day + 1)

    def sum_trials(self,data, group, keep, key):
        '''Get total number of trials grouped by given group. Only use on aggregated
        trial data.
        
        Parameters:
        data --- data frame containing trial data. Must contain all columns in
                `group`, `keep`, and `key`.
        group --- list of columns to group by
        keep --- list of columns to keep
        key --- column containing trial counts'''
        tot = data.groupby(group)[key].sum().rename('trial no')
        keep = data.groupby(group).first()[keep]
        res = pd.concat([tot, keep], axis=1).reset_index()
        return res

    def align_to_timebin(self,data, time_bin, key_dict, keep, index, values):
        '''Align and bin trial data to timebin. Labels each sample with a 
        timebin as a timedelta and float, but does not perform any
        aggregation.
        
        Parameters:
        data --- data frame containing trial data.
        time_bin --- length of bin in minutes to group data by
        key_dict --- names of relevant columns in data
        keep --- list of columns to keep
        index --- list of columns to group by
        values --- list of columns that will be aggregated
        '''
        time = key_dict["timestamp"]
        acc = key_dict["acc"]
        millisbin = key_dict['millis bin']
        offset = key_dict['offset']

        key_start = "trial start"
        key_delta = "delta"

        # label samples with timebin as timestamp
        data = self.get_first_sample(data, index, time, key_start)
        if millisbin != None:
            index.append(millisbin)
            keep.remove(millisbin)
        data = self.bin_by_time(data, time_bin, "hr", index, [], keep + [time, offset] + values, key_start,
                        offset=pd.to_timedelta(12, unit="h"))
        if millisbin != None:
            index.pop()
            keep.append(millisbin)

        # calculate timebin from start of sensory training
        data = self.delta(data, index, key_start, acc, offset, key_delta)

        data = self.time_to_float(data, "Time (hr)", key_delta, "m")
        data["Time (hr)"] = data["Time (hr)"]/60.
        return data

    def align_to_day(self,data, key_dict, keep, index, values):
        '''Align and bin trial data to day. Labels each sample with a 
        day as a timedelta and float, but does not perform any
        aggregation.
        
        Parameters:
        data --- data frame containing trial data.
        time_bin --- length of bin in minutes to group data by
        key_dict --- names of relevant columns in data
        keep --- list of columns to keep
        index --- list of columns to group by
        values --- list of columns that will be aggregated
        '''
        time = key_dict["timestamp"]
        acc = key_dict["acc"]
        millisbin = key_dict['millis bin']
        offset = key_dict['offset']

        key_start = "trial start"
        key_delta = "delta"
        key_day = "day_delta"

        # label samples with day as timestamp
        data = self.get_first_sample(data, index, time, key_start)
        if millisbin != None:
            index.append(millisbin)
            keep.remove(millisbin)
        data = self.bin_by_time(data, 1, "day", index, [], 
                        keep + values + [time, offset, key_delta, "Time (hr)"], key_start, 
                        offset=pd.to_timedelta(12, unit="h"))
        if millisbin != None:
            index.pop()
            keep.append(millisbin)

        # calculate day from start of sensory training
        data = self.delta(data, index, key_start, acc, offset, key_day)

        data = self.time_to_float(data, "Day", "day_delta", "D")
        
        data = data.drop(columns=key_start)
        return data

    def drop_group(self,data, mins, key_dict, key, index, keep):
        '''Drops bins with fewer than the given number of trials and returns data 
        without those bins.

        Filters bins based on total number of samples in index-based groups, 
        as well as number of samples from given key in a group. This key must be 
        binary (e.g. stimulus or blank; water or no water) for the key specific 
        functionality to work. That is, cond0 and cond1 must be from the same
        column, mutually exclusive, and have no other possible values.

        Parameters:
        data --- data frame containing trial data. Must contain columns in `index`,
                `col`, and `key`.
        mins --- dictionary containing
            `min_trials` --- minimum total number of trials to keep a bin (exclusive)
            `min_stimulus` --- minimum number of stimulus trials
            `min_blank` --- minimum number of blank trials
        key --- name of column containing trial number
        index --- list of columns to group data by
        keep --- list of columns to keep
        '''
        kp = keep.copy()
        # include time bin for aggregate grouping    
        index_timebin = index + [key_dict["time bin"]]
        if key_dict["time bin"] in kp:
            kp.remove(key_dict["time bin"])

        # drop bins with fewer than the given number of total, water, or blank trials
        stimulus = key_dict["stimulus type"]
        cond0, cond1 = np.sort(data[stimulus].unique())
        group = data.groupby(index_timebin + [stimulus])
        data = data.set_index(index_timebin)

        #filter bins with fewer than min_blank blank or min_water water trials
        counts = group[key].nunique().unstack(level=stimulus).fillna(0)
        cond = (counts.loc[:, cond0] >= mins["min_blank"]) & (counts.loc[:, cond1] >= mins["min_stimulus"]) 
        data = data[cond]

        #filter bins with fewer than min_trials total trials
        total_group = data.groupby(index_timebin) 
        total_counts = total_group[key].nunique()
        data = data[total_counts >= mins["min_trials"]]
        data = data.reset_index()
        return data

    def agg(self, data, key_dict, index, keep, values):
        '''Aggregate lick frequency values across timebin and day.
        
        Parameters:
        data --- '''
        kp = keep.copy()
        millisbin = key_dict["millis bin"]
        millis = key_dict["Time (ms)"]
        stimulus = key_dict["stimulus type"]
        water = key_dict["water delivery"]

        index_timebin = index + [key_dict["time bin"]]
        if key_dict["time bin"] in kp:
            kp.remove(key_dict["time bin"])

        # number of trials for each timebin by stimulus type
        index_groups = index_timebin + [stimulus, water]
        cond0, cond1 = np.sort(data[stimulus].unique())
        if stimulus in kp:
            kp.remove(stimulus)
        if water in kp:
            kp.remove(water)
        
        counts_groups = self.trial_counts(data, index_groups, kp, key_dict["trial no"])
        
        # mean licking frequencies
        if millisbin != None:
            index_timebin = index_timebin + [millisbin]
            if millisbin in kp:
                kp.remove(millisbin)
                
        data_mean = self.mean_bin(data, index_timebin + [stimulus], values, kp + [water])

        # performance (Lstim - Lblank)
        perf = self.performance(data_mean, index_timebin, kp, stimulus, key_dict["licks"], 
                        cond0, cond1)
        
        return (data_mean, counts_groups, perf)#, poke_perf)

    ############################
    # Main Format and Analysis #
    ############################
    def make_condition_df_helper(self, i=0, animals=[]):
        if i == len(self.andirs):
            for a in animals:
                if not isinstance(a, pd.DataFrame):
                    self.after(1, (lambda : self.evaluate_format(None)))
            animals = pd.concat(animals, ignore_index=True)
            animals["condition"] = self.cond_name
            self.after(1, (lambda : self.evaluate_format(animals)))
            return
        if i < len(self.andirs):
            animal_name = self.andirs[i]
            animal_path = self.in_fold + "\\" + animal_name
            animal = self.make_animal_df(animal_path, animal_name)
            animals.append(animal)
            self.after(1, lambda : self.make_condition_df_helper(i + 1, animals))

    def make_condition_df(self):
        '''Load and format all data files for all animals in a condition, and return as a data frame.
        
        Requires all animals to be included to be in the same directory, which should have no other folders or files in it.
        
        Parameters:
        condir --- path of directory containing animal subdirectories which contain data files
        condition --- name of condition
        metadata --- data frame containing relevant metadata about animals in condition (can contain other information as well)
        meta_vals --- list of colums to include metadata from. Must include 'acc' for length of acclimation (in days).'''
        self.andirs = os.listdir(self.in_fold)
        if len(self.andirs) == 0:
            print('No animals in condition.')
            self.mssg.set('No animals in condition.')
            return pd.DataFrame()
        animals = []
        self.after(1, lambda : self.make_condition_df_helper(i=0, animals=animals))

    def run_formatter(self):
        '''
        Returns formatted training data as a data frame, and writes that data to the
        output file if provided. Can format one or more condition directiories, or one or more animals.
        If multiple animals are being loaded at the same time and `names` is length 1,
        all animals are given the same condition name. If adding to a file of previously
        formatted data, will take the set of columns in the previously formatted data.
        '''
        self.mssg.set(self.mssg.get() + "\nFormatting")
        self.keep = self.keep + ['acc']
        
        # treat provided folder as condition
        if self.as_cond:
            self.after(1, lambda : self.make_condition_df())

        # treat provided folder as animal
        else:
            an_name = self.in_fold.split('\\')[-1].split('/')[-1]
            self.after(1, lambda : self.make_animal_df(self.in_fold, an_name))
        
    def evaluate_format(self,formatted_res):
        if not isinstance(formatted_res, pd.DataFrame):
            self.mssg.set(self.mssg.get() + "\nFormatting failed!")
            self.config(foreground='red')
            self.after(1, (lambda : root.after(2000, (lambda: self.destroy()))))
        else:
            self.data = formatted_res
            self.after(1, (lambda : self.after_format_wrapper()))

    def after_format_wrapper(self):
        keys = self.full_keys.copy()
        keys['millis bin'] = None
        self.fixed_keys = keys
        self.keep = self.keep + ['stimulus', 'water', 'type'] 

        # ignore columns provided in 'keep' that are not present in data
        for c in self.keep:
            if c not in self.data.columns:
                self.keep.remove(c)  

        self.after(1, lambda : self.inst_analyze_wrapper())

    def inst_analyze_wrapper(self):
        if self.ifinst:
            self.mssg.set(self.mssg.get() + "\nInstentaneous analysis")
            self.inst_start() #freq_window, freq_bin, time_bin
        else:
            self.after(1, lambda : self.fixed_window_wrapper())
    
    def inst_start(self):
        self.analysis = self.data.copy()
        self.after(1, lambda: self.inst_time_to_airpuff())

    def inst_time_to_airpuff(self):
        # get time to airpuff delivery
        data = self.puff_delta(self.analysis, ['condition', 'animal','trial no'], "timestamp", "delay", "puff delta")
        
        # correct for pandas behavior when resampling with negative timedeltas
        # if the minimum delay (800 ms) is not present in the dataset
        l = len(data.index)
        if (data['puff delta'].min()) != pd.to_timedelta('-800 ms'):
            data.loc[l, 'trial no'] = 0
            data.loc[l, 'puff delta'] = pd.to_timedelta('-800ms')
        self.analysis = data
        self.after(1, lambda: self.inst_get_lickfreq())

    def inst_get_lickfreq(self):
        # calculate instentanteous or rolling window lick frequency for each trial
        self.analysis = self.bin_by_time(self.analysis, 100, "ms", self.index_trial, self.values, self.keep + ["timestamp", "offset"], "puff delta", origin=None, fn="lickfreq")
        # correct for pandas behavior when resampling with negative timedeltas
        # if the minimum delay (800 ms) is not present in the dataset
        self.analysis = self.analysis.loc[self.analysis['trial no'] > 0]
        self.after(1, lambda: self.inst_align_sample())

    def inst_align_sample(self):
        # align sample times for cross-trial averaging
        self.analysis = self.bin_by_time(self.analysis, 100, "ms", self.index_trial, [], self.keep + ["timestamp", "offset"] + self.values, "puff delta", origin=None)
        self.analysis = self.time_to_float(self.analysis, "Time (ms)", "puff delta", "ms")
        self.after(1, lambda: self.inst_group_timebin())
        
    def inst_group_timebin(self):
        # align to timebin and day
        self.analysis = self.align_to_timebin(self.analysis, self.time_bin, self.full_keys, self.keep + ["puff delta","Time (ms)"], self.index_trial, self.values)
        self.analysis = self.align_to_day(self.analysis, self.full_keys, self.keep + ["puff delta","Time (ms)"], self.index_trial, self.values)
        self.after(1, lambda: self.aggregate_values("full"))

    def fixed_window_wrapper(self):
        if self.iffixed:
            self.mssg.set("Fixed window analysis")
            self.fixed_analysis = self.fixed_window_start() #r, time_bin
        else:
            print("here2")
            self.after(1, lambda: self.finish_analysis())
    
    def fixed_window_start(self):
        self.after(1, lambda: self.fixed_puff_delta())

    def fixed_puff_delta(self):        
        self.fixed_analysis = self.data.copy()
        self.fixed_analysis = self.puff_delta(self.fixed_analysis, ['condition', 'animal', 'trial no'], "timestamp", "delay", "puff delta")
        self.after(1, lambda: self.fixed_first_sample())

    def fixed_first_sample(self):
        # using first sample of each trial to align to timebin and day aligns trial
        # to bin where trial started in edge case where trial spans bins and 
        # anticipatory period happens after bin end
        # matches Matlab v16 analysis behavior
        self.fixed_analysis = self.get_first_sample(self.fixed_analysis, ['condition', 'animal', 'trial no'], 'timestamp', "first sample")
        self.fixed_keys['timestamp'] = 'first sample'
        root.after(1, lambda: self.fixed_lick_freq())

    def fixed_lick_freq(self):
        # get lick frequency for each trial in given window
        self.fixed_analysis = self.fixed_window_get_lickfreq(self.fixed_analysis, ['condition', 'animal', 'trial no'], self.values, 
                                self.keep + ['timestamp', 'first sample','offset'],'puff delta', self.fixedstart, self.fixedend)
        self.after(1, lambda: self.fixed_group_timebin())

    def fixed_group_timebin(self):
        # align trials to timebin and day
        self.fixed_analysis = self.align_to_timebin(self.fixed_analysis, self.time_bin, self.fixed_keys, self.keep, self.index_trial, self.values)
        self.fixed_analysis = self.align_to_day(self.fixed_analysis, self.fixed_keys, self.keep, self.index_trial, self.values)
        self.fixed_keys['timestamp'] = "timestamp"
        self.after(1, lambda: self.aggregate_values("bin"))

    def nth_part(self):
        if self.ifnthfrac:
            self.mssg.set("Last 20% analysis")
            self.nth_part_anal = self.get_nth_part_day(self.fixed_analysis, self.denom, self.num) #num_parts, nth_part
            self.aggregate_values("nth_part")
        else:
            self.after(1, lambda: self.finish_analysis())

    def write_analysis(self,out_files, new_data, columns, add=False):
        mode = 'a' if add else 'w'
        
        for k in range(len(new_data)):
            if add:
                with open(out_files[k]) as f:
                    old_cols = f.readline().strip('\n').split(',')
                out = new_data[k][old_cols]
                cols = old_cols
            else:
                out = new_data[k][columns[k]]
                cols = columns[k]
            out.to_csv(out_files[k], index=False, mode=mode, header=(not add), columns=cols)
            
    def aggregate_values(self, kind):
        all_cols = self.index + self.keep + ["Day", "Time (hr)"]
        if kind =='full':
            fs = [n for n in self.out_files if "inst" in n]
            all_cols = all_cols + ["Time (ms)"]
            keys = self.full_keys
            data = self.analysis
            keys["time bin"] = "delta"
        elif kind == 'bin':
            fs = [n for n in self.out_files if "fixed" in n]
            keys = self.fixed_keys
            data = self.fixed_analysis
            keys["time bin"] = "delta"
        elif kind == 'nth_part':
            fs = [n for n in self.out_files if "part" in n]
            keys = self.fixed_keys
            data = self.nth_part_anal
            keys["time bin"] = "day_delta"
        fs.sort()
        
        keep = list(data.columns) 
        print(keep)
        for c in self.index + self.values + ['trial no']:
            if c in keep:
                keep.remove(c)
        means, counts, perf = self.agg(data, keys, self.index, keep, self.values)
        
        data_cols = all_cols + ['trial no', 'delta', 'day_delta'] + self.values 
        means_cols = all_cols + self.values
        perf_cols = means_cols.copy()
        for i in ['stimulus', 'water', 'type']:
            if i in perf_cols:
                perf_cols.remove(i)
        trial_cols = all_cols + ['trial no']
        if kind =='nth_part':
            if fs:
                columns = [data_cols, means_cols, perf_cols, trial_cols]
                self.write_analysis(fs,[data, means, perf, counts], columns, add=self.add)
            self.after(1, lambda: self.finish_analysis())
            return
        
        daykeep = list(counts.columns)
        for c in self.index + ['Day', 'trial no', 'timestamp', 'offset']:
            if c in daykeep:
                daykeep.remove(c)

        daytots = self.sum_trials(counts, self.index + ["Day"], daykeep, "trial no")
        data_filtered = self.drop_group(data, self.mins, keys, "trial no", self.index, keep)
        filtmeans, filtcounts, filtperf = self.agg(data_filtered, keys, self.index, keep, self.values)
        
        if kind == 'full':
            columns = [data_cols,means_cols, means_cols, perf_cols, perf_cols]        
            output = (data, means, filtmeans, perf,  filtperf)
        elif kind=='bin':
            columns = [data_cols,means_cols, means_cols, perf_cols, perf_cols, trial_cols, trial_cols]        
            output = (data, means, filtmeans, perf,  filtperf, counts, daytots)
        if fs: 
            self.write_analysis(fs,output, columns, add=self.add)
        if kind == "full":
            self.after(1, lambda : self.fixed_window_wrapper())
        elif kind=='bin':
            self.after(1, lambda: self.nth_part())
                 
    def close_and_plot(self):
        if self.plot:
            self.after(1, lambda: self.after(5000, self.destroy))
            self.after(1, lambda: self.master.generate_plots())
        else:
            self.after(1, lambda: self.after(5000, self.destroy))

    def finish_analysis(self):
        self.mssg.set("Done!")
        self.after(1, lambda: self.close_and_plot())
  
class AnalysisApp(ResizeEqualLabelFrame):
    def __init__(self, master=None,**kwargs):
        master.title("Behavior Analysis") #Makes the title that will appear in the top left
        master.bind_all('<Escape>', lambda e: e.widget.focus_set())
        ResizeEqualLabelFrame.__init__(self, master=master,cols=2, rows=2, **kwargs)
        
        # Widgets for analysis
        self.analysis_lf = ResizeEqualLabelFrame(self, rows=5, cols=1,text="1. Analysis")

        self.files_lf = ResizeEqualLabelFrame(self.analysis_lf, rows=2, cols=8,text='Select files for analysis')
        self.flermsg = StringVar()
        self.flerlab = WrappingLabel(self.files_lf, textvariable=self.flermsg, foreground='red')

        self.files_var = StringVar(value='')
        self.files_butt = EnterButton(self.files_lf, text='Select folder(s)', command=(lambda : self.open_directory(self.files_var, self.files_err)))
        self.files_lab = WrappingLabel(self.files_lf, textvariable=self.files_var)
        #files_lab.bind('<Configure>', lambda e: files_lab.config(width=files_lf.winfo_width()))
        self.files_errmsg = StringVar()
        self.files_err = WrappingLabel(self.files_lf, textvariable=self.files_errmsg, foreground='red')
        
        self.metadata_var = StringVar(value='')
        self.metadata_butt = EnterButton(self.files_lf, text='Select metadata file', command=(lambda : self.open_file(self.metadata_var, self.metadata_err)))
        self.metadata_lab = WrappingLabel(self.files_lf, textvariable=self.metadata_var, width=50)
        self.metadata_errmsg = StringVar()
        self.metadata_err = WrappingLabel(self.files_lf, textvariable=self.metadata_errmsg, foreground='red')

        self.formatimemsg = 'Time bin must be a whole number (min 1)'
        self.errtimemsg = StringVar()
        self.time_bin = StringVar(value='4')
        self.check_whole_time_wrapper = (master.register(self.check_whole_time), '%P')
        self.time_en = ttk.Entry(self.files_lf, textvariable=self.time_bin,validate='focusout', validatecommand=self.check_whole_time_wrapper)
        self.time_lab = WrappingLabel(self.files_lf,text='Time bin (hours):')
        self.time_err = WrappingLabel(self.files_lf, textvariable=self.errtimemsg, foreground='red')\

        self.keep_var = StringVar(value='Age, Cage, Sex')
        self.keep_en = ttk.Entry(self.files_lf, textvariable=self.keep_var)
        self.keep_lab = WrappingLabel(self.files_lf, text='Metadata columns to keep:')

        # TODO: validate values in 'lick', 'poke'
        self.errvalmsg = StringVar()
        self.values_var = StringVar(value='lick')
        self.values_en = ttk.Entry(self.files_lf, textvariable=self.values_var)
        self.values_lab = WrappingLabel(self.files_lf,text='Values to analyze:')
        self.values_err = WrappingLabel(self.files_lf, textvariable=self.errvalmsg, foreground='red')

        self.analy_lf = ResizeEqualLabelFrame(self.analysis_lf,rows=2, cols=10, text='Analysis Types:')

        self.analermsg = StringVar()
        self.analerlab = WrappingLabel(self.analy_lf, textvariable=self.analermsg, foreground='red')
        self.inst_anal = BooleanVar(value=False)
        self.inst_cbutt = EnterCheckbutton(self.analy_lf, text='Instentaneous analysis',variable=self.inst_anal, onvalue=True, offvalue=False)

        self.fixed_anal = BooleanVar(value=True)
        self.fixed_cbutt = EnterCheckbutton(self.analy_lf, text='Fixed window (anticipatory) analysis',
                                            variable=self.fixed_anal, onvalue=True, offvalue=False, command=self.enableFixedWindow)

        self.errmsg = StringVar()
        self.formatmsg = 'Start and end should be numbers between -800 and 2000'
        self.check_start_wrapper = (master.register(self.check_start), '%P')
        self.check_end_wrapper = (master.register(self.check_end), '%P')
        
        self.start_val = StringVar(value='700')
        self.startlabel = WrappingLabel(self.analy_lf, text='Start of fixed window:')
        self.fixed_anal_start_en = ttk.Entry(self.analy_lf,textvariable=self.start_val, validate='focusout', validatecommand=self.check_start_wrapper)

        self.end_val = StringVar(value='1000')
        self.endlabel = WrappingLabel(self.analy_lf, text='End of fixed window:')

        self.fixed_anal_end_en = ttk.Entry(self.analy_lf,textvariable=self.end_val, validate='focusout', validatecommand=self.check_end_wrapper)
        self.fixed_anal_err = WrappingLabel(self.analy_lf, textvariable=self.errmsg, foreground='red')

        self.nth_frac_val = BooleanVar(value=True)
        self.nthfrac_cbutt = EnterCheckbutton(self.analy_lf, text='Fraction of day analysis (default last 20%)', 
                                              variable=self.nth_frac_val, onvalue=True, offvalue=False, command=self.enableNthFrac)
        
        self.nthfrac_errmsg = StringVar()
        self.num_val = StringVar(value='5/5')
        self.check_frac_wrapper = (master.register(self.check_frac), '%P')
        self.en_num = ttk.Entry(self.analy_lf, textvariable=self.num_val, validate='focusout', validatecommand=self.check_frac_wrapper)
        self.nth_frac_lab = WrappingLabel(self.analy_lf, text='Fraction of day:')
        self.nth_frac_err = WrappingLabel(self.analy_lf, textvariable=self.nthfrac_errmsg, foreground='red')

        self.formattrialmsg = 'Trial thresholds must be a whole number.'
        self.errtrialmesg = StringVar(value='')
        self.trial_thresh_lf = ResizeEqualLabelFrame(self.analysis_lf, rows=2, cols=6,text='Minimum trial thresholds')
        self.check_whole_wrapper = (master.register(self.check_whole), '%P')
        self.min_var = StringVar(value='10')
        self.min_en = ttk.Entry(self.trial_thresh_lf, textvariable=self.min_var, validate='focusout', validatecommand=self.check_whole_wrapper)
        self.min_lab = WrappingLabel(self.trial_thresh_lf, text='Minimum number of trials:')

        self.min_blank_var = StringVar(value='1')
        self.min_bl_en = ttk.Entry(self.trial_thresh_lf, textvariable=self.min_blank_var, validate='focusout', validatecommand=self.check_whole_wrapper)
        self.min_bl_lab = WrappingLabel(self.trial_thresh_lf, text='Minimum number of blank trials:')

        self.min_stimulus_var = StringVar(value='1')
        self.min_st_en = ttk.Entry(self.trial_thresh_lf, textvariable=self.min_stimulus_var, validate='focusout', validatecommand=self.check_whole_wrapper)
        self.min_st_lab = WrappingLabel(self.trial_thresh_lf, text='Minimum number of stimulus trials:')
        self.min_err = WrappingLabel(self.trial_thresh_lf, textvariable=self.errtrialmesg, foreground='red')


        self.add_lf = ResizeEqualLabelFrame(self.analysis_lf, cols=2, rows=3,text='Output folder')

        self.add_file_label_var = StringVar(value='')
        self.add_file_errmsg = StringVar()
        self.add_file_err = WrappingLabel(self.add_lf, textvariable=self.add_file_errmsg, foreground='red')
        self.add_file_butt = EnterButton(self.add_lf, text='Select output folder',command=(lambda : self.open_directory(self.add_file_label_var, self.add_file_err)))
        self.add_file_label = WrappingLabel(self.add_lf, textvariable=self.add_file_label_var, width=50)

        self.add_val = BooleanVar(value=False)
        self.add_cbutt = EnterCheckbutton(self.add_lf, text='Add to previous analysis',variable=self.add_val, onvalue=True, offvalue=False)

        self.type_lf = ResizeEqualLabelFrame(self.analysis_lf,rows=2, cols=4, text='File type')

        self.type_var = StringVar(value='Animal')
        self.animal_rad = ttk.Radiobutton(self.type_lf, text='Animal', variable=self.type_var, value='Animal', command=self.updateName)
        self.condition_rad = ttk.Radiobutton(self.type_lf, text='Condition', variable=self.type_var, value='Condition',command=self.updateName)

        self.name_lab_var = StringVar()
        self.updateName()
        self.name_lab = WrappingLabel(self.type_lf, textvariable=self.name_lab_var)
        self.name_var = StringVar()
        self.name_en = ttk.Entry(self.type_lf, textvariable=self.name_var)

        self.run_anal_butt = ttk.Button(self.analysis_lf, text='Analyze', command=self.analysis)

        self.plot_lf = ResizeEqualLabelFrame(self,cols=1, rows=5, text="2. Plotting")

        self.load_lf = ResizeEqualLabelFrame(self.plot_lf,cols=2, rows=6, text='Data to plot')
        self.loadfilevar = StringVar()
        self.loadfileerrmsg = StringVar()
        self.loadfileerr = WrappingLabel(self.load_lf, textvariable=self.loadfileerrmsg, foreground='red')
        self.loadbutt = EnterButton(self.load_lf, text='Select Folder', command=(lambda :self.open_directory(self.loadfilevar, self.loadfileerr)))
        self.loadfile = WrappingLabel(self.load_lf,textvariable=self.loadfilevar)

        self.ant_lf = ResizeEqualLabelFrame(self.plot_lf, rows=2, cols=1,text='Fixed Window Plots')
        self.ant_lf_var = BooleanVar(value=True)
        self.ant_lf_cb = EnterCheckbutton(self.ant_lf, text='Fixed window lick frequency', variable=self.ant_lf_var,onvalue=True, offvalue=False)
        self.ant_perf_var = BooleanVar(value=True)
        self.ant_perf = EnterCheckbutton(self.ant_lf, text='Fixed window performance', variable=self.ant_perf_var,onvalue=True, offvalue=False)

       # self.analplot = BooleanVar(value=False)
        
        self.full_lf = ResizeEqualLabelFrame(self.plot_lf,rows=4, cols=1, text='Instentaneous Plots')

        self.fulllckhr = BooleanVar(value=True)
        self.fulllckhrcb = EnterCheckbutton(self.full_lf, text='Instentaneous lick frequency by hour', variable=self.fulllckhr,onvalue=True, offvalue=False)
        self.fullpfhr = BooleanVar(value=True)
        self.fullpfhrcb = EnterCheckbutton(self.full_lf, text='Instentaneous performance by hour', variable=self.fullpfhr,onvalue=True, offvalue=False)
        self.htmplck = BooleanVar(value=True)
        self.htmplckcb = EnterCheckbutton(self.full_lf, text='Instentaneous lick frequency heatmap', variable=self.htmplck,onvalue=True, offvalue=False)
        self.htmppf = BooleanVar(value=True)
        self.htmppfcb = EnterCheckbutton(self.full_lf, text='Instentaneous performance heatmap', variable=self.htmppf,onvalue=True, offvalue=False)

        self.trial_lf = ResizeEqualLabelFrame(self.plot_lf,rows=3, cols=1, text='Trial Plots')
        self.last20val = BooleanVar(value=True)
        self.last20cb = EnterCheckbutton(self.trial_lf, text='Nth fraction of day stimulus vs blank', variable=self.last20val,onvalue=True, offvalue=False)
        self.trialsbyhour = BooleanVar(value=True)
        self.trialsbyhourcb = EnterCheckbutton(self.trial_lf, text='Trials by hour', variable=self.trialsbyhour, onvalue=True, offvalue=False)

        self.other_lf = ResizeEqualLabelFrame(self.plot_lf,cols=1,rows=1, text='Other Plots')

        self.save_lf = ResizeEqualLabelFrame(self.plot_lf,cols=2, rows=4, text='Save Options')
        self.savevar = BooleanVar(value=True)
        self.savecb = EnterCheckbutton(self.save_lf, text='Save plots', variable=self.savevar,onvalue=True, offvalue=False, command=(lambda:self.disablesavebutt()))

        self.savefilevar = StringVar()
        self.savefileerrmsg = StringVar()
        self.savefileerr = WrappingLabel(self.save_lf, textvariable=self.savefileerrmsg, foreground='red')
        self.savebutt = EnterButton(self.save_lf, text='Select Folder', command=(lambda : self.open_directory(self.savefilevar, self.savefileerr)))
        self.savelab = WrappingLabel(self.save_lf, textvariable=self.savefilevar, width=25)
    
        self.dispvar = BooleanVar(value=True)
        self.dispcb = EnterCheckbutton(self.save_lf, text='Display plots', variable=self.dispvar,onvalue=True, offvalue=False)

        self.run_plots_butt = EnterButton(self.plot_lf, text='Generate plots', command=self.plot)
        self.analyplot = BooleanVar(value=False)
        self.run_all_butt = EnterButton(self, text='Analyze and plot', command=self.analyze_and_plot)

        self.defaultpad=2
        self.defoutpad=5
        self.grid_all()

    def grid_all(self):
                # grid plot buttons
        self.plot_lf.grid(row=0,column=1, padx=self.defoutpad, pady=self.defoutpad,sticky='nsew')

        self.load_lf.grid(row=0,column=0, padx=self.defoutpad, pady=self.defoutpad,sticky='nsew')
        self.loadbutt.grid(row=0,column=0, padx=self.defaultpad, pady=self.defaultpad,sticky='nsew')
        self.loadfile.grid(row=0,column=1, padx=self.defaultpad, pady=self.defaultpad,sticky='nsew')

        self.ant_lf.grid(row=1,column=0, padx=self.defoutpad, pady=self.defoutpad,sticky='nsew')
        self.ant_perf.grid(row=1,column=0, padx=self.defaultpad, pady=self.defaultpad,sticky='nsew')
        self.ant_lf_cb.grid(row=0,column=0, padx=self.defaultpad, pady=self.defaultpad,sticky='nsew')

        self.full_lf.grid(row=2,column=0, padx=self.defoutpad, pady=self.defoutpad,sticky='nsew')
        self.fullpfhrcb.grid(row=1,column=0, padx=self.defaultpad, pady=self.defaultpad,sticky='nsew')
        self.fulllckhrcb.grid(row=0,column=0, padx=self.defaultpad, pady=self.defaultpad,sticky='nsew')
        self.htmppfcb.grid(row=3,column=0, padx=self.defaultpad, pady=self.defaultpad,sticky='nsew')
        self.htmplckcb.grid(row=2,column=0, padx=self.defaultpad, pady=self.defaultpad,sticky='nsew')

        self.trial_lf.grid(row=3,column=0, padx=self.defoutpad, pady=self.defoutpad,sticky='nsew')
        self.last20cb.grid(row=0,column=0, padx=self.defaultpad, pady=self.defaultpad,sticky='nsw')
        self.trialsbyhourcb.grid(row=1,column=0, padx=self.defaultpad, pady=self.defaultpad,sticky='nsew')

        #self.other_lf.grid()

        self.save_lf.grid(row=6,column=0, padx=self.defoutpad, pady=self.defoutpad,sticky='nsew')
        self.savebutt.grid(row=1,column=0, padx=self.defaultpad, pady=self.defaultpad,sticky='nsew')
        self.savelab.grid(row=1,column=1, padx=self.defaultpad, pady=self.defaultpad,sticky='nsew')
        self.savecb.grid(row=0,column=0, padx=self.defaultpad, pady=self.defaultpad,sticky='nsew')
        self.dispcb.grid(row=3,column=0, padx=self.defaultpad, pady=self.defaultpad,sticky='nsew')

        self.run_plots_butt.grid(row=7,column=0, padx=self.defoutpad, pady=self.defoutpad,sticky='nsew')

        # grid buttons for analysis column
        self.analysis_lf.grid(row=0, column=0, padx=self.defoutpad, pady=self.defoutpad,sticky='ewns')
        self.files_lf.grid(row=0, column=0, padx=self.defaultpad, pady=self.defaultpad,sticky='ewns')

        self.files_butt.grid(row=0, column=0, sticky='nsew',padx=self.defaultpad, pady=self.defaultpad)
        self.files_lab.grid(row=0, column=1, sticky='nsew', padx=self.defaultpad, pady=self.defaultpad)

        self.metadata_butt.grid(row=2, column=0, sticky='nsew',padx=self.defaultpad, pady=self.defaultpad)
        self.metadata_lab.grid(row=2, column=1, sticky='nsew',padx=self.defaultpad, pady=self.defaultpad)

        self.keep_en.grid(row=7, column=1, sticky='nsew',padx=self.defaultpad, pady=self.defaultpad)
        self.keep_lab.grid(row=7, column=0,sticky='ewns',padx=self.defaultpad, pady=self.defaultpad)

        self.values_en.grid(row=8, column=1, sticky='nsew',padx=self.defaultpad, pady=self.defaultpad)
        self.values_lab.grid(row=8, column=0, sticky='ewns',padx=self.defaultpad, pady=self.defaultpad)

        self.time_en.grid(row=5, column=1, padx=self.defoutpad, pady=self.defoutpad,sticky='nsew')
        self.time_lab.grid(row=5, column=0, padx=self.defoutpad, pady=self.defoutpad,sticky='ewns')

        self.analy_lf.grid(row=3, column=0, padx=self.defoutpad, pady=self.defoutpad,sticky='nsew')

        self.inst_cbutt.grid(row=1, column=0,sticky='nsew',padx=self.defaultpad, pady=self.defaultpad, columnspan=2)

        self.fixed_cbutt.grid(row=2, column=0, sticky='nsew',padx=self.defaultpad,pady=self.defaultpad, columnspan=2)
        self.startlabel.grid(row=3, column=0,padx=[50,self.defaultpad], pady=self.defaultpad,sticky='ewns')
        self.fixed_anal_start_en.grid(row=3, column=1,padx=self.defaultpad, pady=self.defaultpad,sticky='w')
        self.endlabel.grid(row=4, column=0,padx=[50,self.defaultpad], pady=self.defaultpad,sticky='nsew')
        self.fixed_anal_end_en.grid(row=4, column=1,padx=self.defaultpad, pady=self.defaultpad,sticky='nsw')

        self.nthfrac_cbutt.grid(row=6, column=0, columnspan=3,padx=self.defaultpad, pady=self.defaultpad,sticky='nsew')
        self.en_num.grid(row=7, column=1,padx=self.defaultpad, pady=self.defaultpad,sticky='w')
        self.nth_frac_lab.grid(row=7, column=0,padx=[50,self.defaultpad], pady=self.defaultpad,sticky='ewns')

        self.trial_thresh_lf.grid(row=4, column=0,padx=self.defoutpad, pady=self.defoutpad,sticky='w')
        self.min_en.grid(row=0, column=1,padx=self.defaultpad, pady=self.defaultpad,sticky='w')
        self.min_lab.grid(row=0, column=0,padx=self.defaultpad, pady=self.defaultpad,sticky='ewns')
        self.min_bl_en.grid(row=1, column=1,padx=self.defaultpad, pady=self.defaultpad,sticky='w')
        self.min_bl_lab.grid(row=1, column=0,padx=self.defaultpad, pady=self.defaultpad,sticky='ewns')
        self.min_st_en.grid(row=2, column=1,padx=self.defaultpad, pady=self.defaultpad,sticky='w')
        self.min_st_lab.grid(row=2, column=0,padx=self.defaultpad, pady=self.defaultpad,sticky='ewns')

        self.add_lf.grid(row=5, column=0,padx=self.defoutpad, pady=self.defoutpad,sticky='w')
        self.add_cbutt.grid(row=2, column=0, columnspan=2,padx=self.defaultpad, pady=self.defaultpad,sticky='ew')
        self.add_file_butt.grid(row=0, column=0,padx=self.defaultpad, pady=self.defaultpad,sticky='nsew')
        self.add_file_label.grid(row=0, column=1,padx=self.defaultpad, pady=self.defaultpad,sticky='nsew')

        self.type_lf.grid(row=6, column=0,padx=self.defaultpad, pady=self.defaultpad, sticky='w')
        self.animal_rad.grid(row=0, column=0, padx=[50,self.defaultpad], pady=self.defaultpad,sticky='nsew')
        self.condition_rad.grid(row=1, column=0,padx=[50,self.defaultpad], pady=self.defaultpad,sticky='nsew')
        self.name_lab.grid(row=2, column=0, padx=[50,self.defaultpad], pady=self.defaultpad,sticky='nse')
        self.name_en.grid(row=2, column=1, padx=self.defaultpad, pady=self.defaultpad,sticky='nsew')

        self.run_anal_butt.grid(row=20, column=0,padx=self.defoutpad, pady=self.defoutpad)

        self.run_all_butt.grid(row=1, column=0, padx=self.defoutpad, pady=self.defoutpad)

    def disablesavebutt(self):
        if self.savevar.get():
            self.savebutt.state(['!disabled'])
        else:
            self.savebutt.state(['disabled'])

    def fixed_num_valid(self, val):
        return re.match('^-?[0-9]*$', val) is not None and int(val) > -800 and int(val) < 2000
    
    def check_end(self, newval):
        valid = self.fixed_num_valid(newval) and self.fixed_num_valid(self.start_val.get()) and int(newval) > int(self.start_val.get())
        if not valid:
            self.fixed_anal_err.grid(row=5, column=0,padx=self.defaultpad, pady=self.defaultpad,sticky='nsw')
            self.errmsg.set(self.formatmsg)
        else:
            self.fixed_anal_err.grid_forget()
        return valid
    
    def check_start(self, newval):
        valid = self.fixed_num_valid(newval) and self.fixed_num_valid(self.end_val.get()) and int(newval) < int(self.end_val.get())
        if not valid:
            self.fixed_anal_err.grid(row=5, column=0,padx=self.defaultpad, pady=self.defaultpad,sticky='nsw')
            self.errmsg.set(self.formatmsg)
        else:
            self.fixed_anal_err.grid_forget()
        return valid

    def check_whole_time(self,newval):
        valid = re.match('^[0-9]*$', newval) is not None and int(newval) > 0
        if not valid:
            self.time_err.grid(row=6, column=0,padx=self.defoutpad, pady=self.defoutpad,sticky='nsw')
            self.errtimemsg.set(self.formatimemsg)
        else:
            self.time_err.grid_forget()
        return valid

    def check_frac(self,newval):
        valid = re.match('^[0-9]/[0-9]$', newval) is not None
        if not valid:
            self.nth_frac_err.grid(row=9, column=0,padx=self.defaultpad, pady=self.defaultpad,sticky='nsw')
            self.nthfrac_errmsg.set('Numerator and denomenator must be digits.')
            return valid
        num = newval.split('/')[0]
        denom = newval.split('/')[1]
        valid = valid and int(num) <= int(denom)
        if not valid:
            self.nth_frac_err.grid(row=9, column=0,padx=self.defaultpad, pady=self.defaultpad,sticky='nsw')
            self.nthfrac_errmsg.set('Numerator must be less than denomenator.')
        else:
            self.nth_frac_err.grid_forget()
        return valid
    
    def check_whole(self,newval):
        valid = re.match('^[0-9]*$', newval) is not None and int(newval) > 0
        if not valid:
            self.errtrialmesg.set(self.formattrialmsg)
            self.min_err.grid(row=3, column=0,padx=self.defaultpad, pady=self.defaultpad,sticky='nsw')
        else: 
            self.min_err.grid_forget()
        return valid

    def updateName(self):
        self.name_lab_var.set(f'{self.type_var.get()} name (optional):')

    def validate_previous_analysis(self):
        self.analermsg.set('')
        self.flermsg.set('')
        self.analerlab.grid_forget()
        self.flerlab.grid_forget()
        prev_anal_files = os.listdir(self.add_file_label_var.get())
        prev_anal_files.sort()
            
        if len(prev_anal_files) == 0:
            self.add_file_errmsg.set("No previous analysis.")
            self.add_file_err.grid(row=1, column=0, padx=self.defaultpad, pady=self.defaultpad, sticky='nsw')
            return
        
        fixed = [n for n in prev_anal_files if 'fixed' in n]
        nth = [n for n in prev_anal_files if 'part' in n]
        full = [n for n in prev_anal_files if 'inst' in n]

        fixed_str = 'fixed, ' if self.fixed_anal.get() else ''
        inst_str = 'inst, ' if self.inst_anal.get() else ''
        prev_inst_str = 'inst, ' if len(full) > 0 else ''
        prev_fixed_str = 'fixed, ' if len(fixed) > 0 else ''
        nth_str = 'part' if len(fixed) > 0 else ''
        prev_nth_str = 'part' if len(nth) > 0 else ''

        if fixed_str != prev_fixed_str or inst_str != prev_inst_str or nth_str != prev_nth_str:
            self.analermsg.set('Current analysis selection (%s %s %s) does not match previous analysis (%s %s %s)' 
                        % (inst_str,fixed_str, nth_str, prev_inst_str, prev_fixed_str, prev_nth_str))
        
        fixed_start, fixed_end = fixed[0].split('_')[0].split(' ')[-1].split('-')

        if fixed_start != self.start_val.get():
            self.analermsg.set(self.analermsg.get() + f'\nCurrent fixed start ({self.start_val.get()}) does not match previous fixed start ({fixed_start})')
        if fixed_end != self.end_val.get():
            self.analermsg.set(self.analermsg.get() + f'\nCurrent fixed end ({self.end_val.get()}) does not match previous fixed end ({fixed_end})')
        
        num, denom = nth[0].split(' ')[0].split('-')
        cnum, cdenom = self.num_val.get().split('/')
        if num != cnum or denom != cdenom:
            self.analermsg.set(self.analermsg.get() +f'\nCurrent frac day ({cnum}/{cdenom}) does not match previous frac day ({num}/{denom})')
        
        testfile = pd.read_csv(self.add_file_label_var.get() + '\\' + prev_anal_files[-1])

        prevcols = list(testfile.columns.drop(['condition', 'animal', 'acc', 'stimulus', 'water', 'type', 'Day', 'Time (hr)', 'Time (ms)', 'trial no', 'lick', 'poke','delta', 'day_delta'], errors='ignore'))
        currcols = self.keep_var.get().split(', ')
        if Counter(prevcols) != Counter(currcols):
            self.flermsg.set(f"Current columns to keep ({', '.join(currcols)}) do not match previous columns to keep ({', '.join(prevcols)})")

        prevtim = testfile.groupby(['condition', 'animal', 'Time (hr)']).first().reset_index()['Time (hr)'].diff(-1).abs()[0]
        currtim = float(self.time_bin.get())
        if not math.isclose(prevtim, currtim):
            self.flermsg.set(self.flermsg.get() + f'\nCurrent time bin ({currtim}) does not match previous timebin ({prevtim})')

        prevval = list(testfile.columns.drop(['condition', 'animal', 'acc', 'stimulus', 'water', 'type', 'Day', 'Time (hr)', 'trial no', 'Time (ms)', 'delta', 'day_delta'] + prevcols, errors='ignore'))
        currval = self.values_var.get().split(', ')
        if Counter(prevval) != Counter(currval):
            self.flermsg.set(self.flermsg.get() + f"\nCurrent values ({', '.join(currval)}) do not match previous values ({', '.join(prevval)})")

        if self.flermsg.get() != '':
            self.flerlab.grid(row=0, column=0, padx=self.defaultpad, pady=self.defaultpad, sticky='nsew')
        if self.analermsg.get() != '':
            self.analerlab.grid(row=0, column=0, padx=self.defaultpad, pady=self.defaultpad, sticky='nsew')
            print(self.analermsg.get())
        if self.flermsg.get() == '' and self.analermsg.get() == '':
            return True
        else:
            return False

    # Function to open a file dialog
    def open_file(self, file_path_var, file_err):
        file_path = filedialog.askopenfilename()
        file_path_var.set(file_path)
        file_err.grid_forget()

    def open_directory(self, file_path_var, file_err):
        file_path = filedialog.askdirectory()
        file_path_var.set(file_path)
        file_err.grid_forget()
        
    def enableFixedWindow(self):
        if self.fixed_anal.get():
            self.fixed_anal_start_en.state(['!disabled'])
            self.fixed_anal_end_en.state(['!disabled'])
            self.nthfrac_cbutt.state(['!disabled'])
            self.en_num.state(['!disabled'])
        else:
            self.fixed_anal_start_en.state(['disabled'])
            self.fixed_anal_end_en.state(['disabled'])
            self.nthfrac_cbutt.state(['disabled'])
            self.en_num.state(['disabled'])

    def enableAddFolder(self):
        if self.add_val.get():
            self.add_file_butt.state(['!disabled'])
        else:
            self.add_file_butt.state(['disabled'])

    def enableNthFrac(self):
        if self.nth_frac_val.get():
            self.en_num.state(['!disabled'])
        else:
            self.en_num.state(['disabled'])

    def analysis(self):
        if self.files_var.get() == '':
            self.files_errmsg.set("Please select data to analyze.")
            self.files_err.grid(row=1, column=0, padx=self.defaultpad, pady=self.defaultpad, sticky='nsw')
            return
        if self.metadata_var.get() == '':
            self.metadata_errmsg.set("Please select metadata file.")
            self.metadata_err.grid(row=3, column=0, padx=self.defaultpad, pady=self.defaultpad, sticky='nsw')
            return
        if self.add_file_label_var.get() == '':
            self.add_file_errmsg.set('Please select output file.')
            self.add_file_err.grid(row=1, column=0, padx=self.defaultpad, pady=self.defaultpad, sticky='nsw')
            return
        if self.add_val.get():
            if not self.validate_previous_analysis():
                return
        if not self.add_val.get() and len(os.listdir(self.add_file_label_var.get())) > 0:
            res = tk.messagebox.askyesno(icon=tk.messagebox.WARNING, message='Analyzing without adding will overwrite any previous analysis in this folder. Proceed?')
            if not res:
                return
        
        RunAnalysisWindow(self)
      
    def plot(self):
        f = self.loadfilevar.get()
        if f == '':
            self.loadfileerrmsg.set("Please select data to plot.")
            self.loadfileerr.grid(row=1,column=1)
            return
        if self.savevar.get() and self.savefilevar.get() == '':
            self.savefileerrmsg.set("Please select a folder to save plots.")
            self.savefileerr.grid(row=2, column=0)
            return
        
        PlotWindow(self)
        
    def analyze_and_plot(self):
        if self.savevar.get() and self.savefilevar.get() == '':
            self.savefileerrmsg.set("Please select a folder to save plots.")
            self.savefileerr.grid(row=2, column=0)
            return
        self.analyplot.set(True)
        self.loadfilevar.set(self.add_file_label_var.get())
        self.analysis()
        return

if __name__ == '__main__':
    # Start the main loop
    root = tk.Tk() #Makes the window
    app = AnalysisApp(root)
    app.grid(row=0, column=0, sticky='nsew')
    root.mainloop()
