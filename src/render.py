from joypy import joyplot
from matplotlib.pyplot import show
from seaborn import catplot
from pandas import DataFrame


def renderCount(dict1: dict, dict2: dict, header1: str, header2: str):
    dataframe = DataFrame.from_records(data=[dict1, dict2])
    dataframe.insert(loc=0, column="Species", value=[header1, header2])
    dataframe = dataframe.melt(id_vars=['Species'], var_name="Amino Acid", value_name="Residues")
    catplot(kind="bar", data=dataframe, hue='Species', y="Residues", x="Amino Acid")
    show()
