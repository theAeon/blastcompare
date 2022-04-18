from typing import Tuple
from joypy import joyplot
from matplotlib.pyplot import show
from seaborn import catplot, pointplot, countplot, barplot, heatmap
from pandas import DataFrame, concat
from matplotlib import pyplot, axes, figure
from mpl_toolkits import mplot3d


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

def rendermulti(renderTuple: list[tuple[dict, str]]):
    figtmp, axitmp = pyplot.subplots(figsize=(9,6))
    combilong = concat([singleFrame(d, s) for d, s in renderTuple])
    combi = combilong.pivot("Species", "Amino Acid", "Residues")
    heatmap(combi, ax=axitmp, yticklabels=False, xticklabels=False)
    
    return figtmp

