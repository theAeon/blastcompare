from typing import Tuple
from matplotlib.pyplot import show
from seaborn import catplot, pointplot, heatmap
from pandas import DataFrame, concat
from matplotlib import pyplot, axes


def renderCount(dict1: dict, dict2: dict, header1: str, header2: str):
    dataframe = DataFrame.from_records(data=[dict1, dict2])
    dataframe.insert(loc=0, column="Species", value=[header1, header2])
    dataframe = dataframe.melt(id_vars=['Species'], var_name="Amino Acid", value_name="Residues")
    catplot(kind="bar", data=dataframe, hue='Species', y="Residues", x="Amino Acid")
    show()

def singleFrame(dict: dict, header: str):
    dataframe = DataFrame.from_records([dict])
    dataframe.insert(loc=0, column="Species", value=[header])
    try:
        dataframe = dataframe.drop(columns=["GAP"])
    except KeyError:
        pass
    except BaseException:
        raise
    dataframe = dataframe.melt(
        id_vars=['Species'], var_name="Amino Acid", value_name="Residues")
    return dataframe

def renderSingle(dict: dict, header: str) -> axes.Axes:
    dataframe = singleFrame(dict, header)
    fig, ax = pyplot.subplots()
    pointplot(data=dataframe, hue='Species',
          y="Residues", x="Amino Acid", ax=ax)
    return ax

class heatmapper():
    renderTuple: list[tuple[dict, str]]
    def __init__(self, renderTup):
        self.dfcol=concat([singleFrame(d, s) for d, s in renderTup])
        self.figtmp, self.axitmp = pyplot.subplots(figsize=(9, 6))
        self.evnum =self.figtmp.canvas.mpl_connect('button_press_event', self.onclick)
    def rendermulti(self):
        combi = self.dfcol.pivot("Species", "Amino Acid", "Residues")
        heatmap(combi, ax=self.axitmp, yticklabels=False)

    def onclick(self, event):
        print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
              ('double' if event.dblclick else 'single', event.button,
               event.x, event.y, event.xdata, event.ydata))
