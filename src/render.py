from math import floor
from matplotlib.pyplot import show
from seaborn import catplot, pointplot, heatmap, light_palette
from pandas import DataFrame, concat
from matplotlib import pyplot, axes
from grid_strategy import strategies


def renderCount(dict1: dict, dict2: dict, header1: str, header2: str):
    dataframe = DataFrame.from_records(data=[dict1, dict2])
    dataframe.insert(loc=0, column="Species", value=[header1, header2])
    dataframe = dataframe.melt(
        id_vars=['Species'], var_name="Amino Acid", value_name="Residues")
    catplot(kind="bar", data=dataframe, hue='Species',
            y="Residues", x="Amino Acid")
    show()


def singleFrame(dicta: dict, header: str):
    dataframe = DataFrame.from_records([dicta])
    dataframe.insert(loc=0, column="Species", value=[header])
    try:
        dataframe = dataframe.drop(columns=["GAP"])
    except KeyError:
        pass
    except BaseException:
        print("fatal")
        raise
    dataframe = dataframe.melt(
        id_vars=['Species'], var_name="Amino Acid", value_name="Residues")
    return dataframe


def renderSingle(dicta: dict, header: str) -> axes.Axes:
    dataframe = singleFrame(dicta, header)
    fig, ax = pyplot.subplots()
    pointplot(data=dataframe, hue='Species',
              y="Residues", x="Amino Acid", ax=ax)
    return ax


class heatmapper():
    renderTuple: list[tuple[dict, str]]
    
    def __init__(self):
        self.combi = None
        self.axitmp = None
        self.figtmp = pyplot.figure(figsize=(17.5, 10.25), layout='tight')
        self.colors = None
        self.evnum = self.figtmp.canvas.mpl_connect(
            'button_press_event', self.onclick)

    def rendermulti(self, cbar=True, vmin=None, vmax=None, center=None, cmap=None):
        heatmap(self.combi, ax=self.axitmp, yticklabels=False, cbar=cbar,
                vmin=vmin, vmax=vmax, center=center,
                cbar_kws={'use_gridspec': True}, cmap=cmap)
        print(self.combi)

    def onclick(self, event):
        print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
              ('double' if event.dblclick else 'single', event.button,
               event.x, event.y, event.xdata, event.ydata))
        print(self.combi.iloc[floor(event.ydata)])

class singleheatmapper(heatmapper):
    def __init__(self, renderTup):
        super().__init__()
        dfcol = concat([singleFrame(d, s) for d, s in renderTup])
        self.combi = dfcol.pivot(
            "Species", "Amino Acid", "Residues")
        self.figtmp, self.axitmp = pyplot.subplots(figsize=(9, 6))



class multiheatmapper(heatmapper):
    def __init__(self, listRenderTup, namelist):
        super().__init__()
        strat = strategies.RectangularStrategy()
        self.spec = strat.get_grid(len(listRenderTup))
        self.listtup = listRenderTup
        self.namelist = namelist
        self.listsub = []
        self.combilist = []
        self.localmin = 0
        self.localmax = 0
        fig = pyplot.figure(1)
        for i, tup in enumerate(self.listtup):
            self.listsub.append(pyplot.subplot(self.spec[i]))
            dftmp = concat([singleFrame(d, s) for d, s in tup])
            self.combilist.append(dftmp.pivot(
            "Species", "Amino Acid", "Residues"))
            num = self.combilist[i].drop(columns=["Val"])
            mi = min(num.min())
            ma = max(num.max())
            self.localmin = mi if mi < self.localmin else self.localmin
            self.localmax = ma if ma > self.localmax else self.localmax
        self.axitmp = None
        self.combi = None
        self.namesel = None
        fig.canvas.mpl_connect('axes_enter_event', self.enter_axes)
        

    def seltmp(self, i):
        self.axitmp = self.listsub[i]
        self.combi = self.combilist[i]
        self.namesel = self.namelist[i]

    def rendermulti(self, cbar=True, vmax=None, vmin=None, ax=None):
        vmax = self.localmax if vmax is None else vmax
        vmin = self.localmin if vmin is None else vmin
        ax = self.axitmp if ax is None else ax
        center = vmin + 85
        cmap = light_palette('red', n_colors=256, as_cmap=True)
        self.axitmp.title.set_text(self.namesel)
        return super().rendermulti(center=center, vmin=vmin, vmax=vmax, cbar=cbar, cmap=cmap)
    
    def enter_axes(self, event):
        print('enter_axes', event.inaxes)
        for i, ax in enumerate(self.listsub):
            if ax.in_axes(event):
                self.seltmp(i)
                print(i)
