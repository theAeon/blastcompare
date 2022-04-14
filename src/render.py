from joypy import joyplot
from matplotlib.pyplot import show
from seaborn import catplot
from pandas import DataFrame


def renderCount(dict1: dict, dict2: dict, header1: str, header2: str):
    dataframe = DataFrame.from_records(data=[dict1, dict2])
    dataframe = dataframe.reset_index(0)
    dataframe = dataframe.replace({0:header1, 1:header2})
    dataframe = dataframe.melt(id_vars=['index'], var_name="Amino Acid", value_name="Length")
    catplot(kind="bar", data=dataframe, hue='index', y="Length", x="Amino Acid")
    show()
