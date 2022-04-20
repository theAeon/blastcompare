from math import floor
from matplotlib.pyplot import show
from seaborn import catplot, pointplot, heatmap
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

    def __init__(self, renderTup):
        self.dfcol = concat([singleFrame(d, s) for d, s in renderTup])
        self.combi = self.dfcol.pivot(
            "Species", "Amino Acid", "Residues")
        self.figtmp, self.axitmp = pyplot.subplots(figsize=(9, 6))
        self.evnum = self.figtmp.canvas.mpl_connect(
            'button_press_event', self.onclick)

    def rendermulti(self):
        heatmap(self.combi, ax=self.axitmp, yticklabels=False)
        print(self.combi)

    def onclick(self, event):
        print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
              ('double' if event.dblclick else 'single', event.button,
               event.x, event.y, event.xdata, event.ydata))
        print(self.combi.iloc[floor(event.ydata)])
        
