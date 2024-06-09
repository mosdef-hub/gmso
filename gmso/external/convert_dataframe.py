"""Module support for converting to/from Pandas DataFrame objects."""

import functools
import warnings
from collections.abc import Iterable
from operator import attrgetter, itemgetter

import numpy as np
import pandas as pd
import unyt as u
from symengine import expand

from gmso import Topology
from gmso.abc.serialization_utils import unyt_to_dict
from gmso.core.views import PotentialFilters

pfilter = PotentialFilters.UNIQUE_PARAMETERS
from gmso.exceptions import GMSOError
from gmso.utils.io import import_

pd = import_("pandas")


def to_dataframeDict(
    topology: Topology,
    parameters: str or list = "all",
    format: str = "default",
    columns: list[str] = None,
    convert_unyts: str = "to_headers",
) -> pd.DataFrame:
    """Return a dictionary of pandas dataframe objects for a topology.

    Parameters
    ----------
    topology : gmso.Topology, required
        Topology to use for converting values
    parameters : str or list of str, optional, default='all'
        A string determining what aspects of the gmso topology will be reported.
        Options are: 'all', 'sites', 'bonds', 'angles', 'dihedrals', and 'impropers'. Defaults to 'all'. Can pass multiple strings as a list.
    format : str, optional, default='default'
        The output formatting style for the dataframe.
        Options are 'default', 'specific_columns', 'publication', `remove_duplicates`. Defaults to 'default'
        'default' will output default column values of ["name", "atom_type.name", "atom_type.parameters", "charge", "mass"],
        and any additional attributes in the `columns` argument.
        'specific_columns' will only output the attributes from the `columns` argument.
        'publication' will use the default outputs, but remove duplicate values from the dataframes. It adds a column labeled
        'Atom Indices' to the `sites` dataframe, which enumerates the indices that the atom_type is a part of.
        `remove_duplicates` will remove duplicate rows from the dataframe. For sites, this column is `atom_types.name`.
        For connections, it is the `connection_types.connection_members`. For sites, an additional column will be added, labeled
        `Atom Indices` that includes the site indexes of members that are the given `atom_type.name`. Because these methods
        are specific to a given Topology element, the `parameters` argument must be one of
        {"sites", "bonds", "angles", "dihedrals", "impropers"}, not {"all"}.
    columns : list of str, optional, default=None
        List of strings that are attributes of the topology site and can be included as entries in the pandas dataframe.
        Examples of these can be found by printing `topology.sites[0].__dict__` or `topology.bonds[0].__dict__`.
        See https://gmso.mosdef.org/en/stable/data_structures.html#gmso.Atom for additional information on labeling.
    convert_unyts: str, optional, default='to_headers'
        The placement/recording of unyt quantities in dataframe.
        Options are 'to_headers', 'with_data', 'no_unyts'
        Determines if numerical values in the DataFrame are saved as unyt quantities or floats. Default case, 'to_headers",
        puts the unyts as strings to go with the column header of the dataframe.
        `with_data` leaves any values alone, so any values in the Topology that are unyt quantities will stay that way.
        `no_unyts` strips any unyt values and converts to a float in the associated element of the dataframe.
        See https://unyt.readthedocs.io/en/stable/usage.html
        for more information about manipulating unyt quantities.

    Returns
    -------
    Dictionary of Pandas Dataframe
        A python dictionary of pandas.Dataframe object, see https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html
        for further information. The keys of this dictionary are the attributes of the topology that are associated with each DataFrame. These
        can be `sites`, `bonds`, `angles`, `dihedrals`, `impropers`, which are determined from the argument `parameters`.


    Examples
    ________
    # example topology to use
    ``` python
    >>> import gmso
    >>> import mbuild as mb
    >>> from gmso.parameterization import apply
    >>> cpd = mb.load("C", smiles=True)
    >>> top = cpd.to_gmso()
    >>> ff = gmso.ForceField("oplsaa")
    >>> ptop = apply(top, ff)
    ```


    >>> gmso.external.convert_dataframe.to_dataframeDict(ptop, parameters='sites', columns=['charge'], convert_unyts="to_headers")
        This will return a dataframe with a listing of the sites and include the charges that correspond to each site.
        ```
        {'sites':
            name atom_type.name  epsilon (kJ/mol)  sigma (nm)   charge (elementary_charge)  mass (amu)
            0    C       opls_138          0.276144        0.35  -0.24      12.011
            1    H       opls_140          0.125520        0.25   0.06       1.008
            2    H       opls_140          0.125520        0.25   0.06       1.008
            3    H       opls_140          0.125520        0.25   0.06       1.008
            4    H       opls_140          0.125520        0.25   0.06       1.008
        }
        ```

    >>> topology.to_dataframe(parameters = 'dihedrals', site_attrs = ['positions'])
        This will return a dataframe with a listing of the sites that make up each dihedral, the positions of each of
        those sites, and the parameters that are associated with the dihedrals.

    Notes
    -----
    A dataframe is easily manipulated. In order to change the rounding to two decimals places for a column named `label`:
        >>> df['label'] = df['label'].round(2)
    The column labels can also be easily modified. This line can take a dataframe `df` and rename a column labeled
    `Atom0` to `newname` using a dictionary.
        >>> df.rename(columns = {'Atom0':'newname'})
    See https://pandas.pydata.org/pandas-docs/stable/getting_started/intro_tutorials/index.html for further information.
    """
    if columns is None:
        columns = []
    if not topology.is_typed():
        raise GMSOError(
            "This topology is not typed, please type this object before converting to a pandas dataframe"
        )
    outDict = {}  # dictionary of dataframes to write out

    # Get columns from format methods
    columnsDict = {}
    connectionsList = [
        "bonds",
        "angles",
        "dihedrals",
        "impropers",
    ]  # these can be handled generally
    remove_duplicatesBool = False  # flag to remove duplicate parameters from the dataframe and put indices into a new column
    if format == "default":
        columnsDict = {
            param: [
                "name",
                f"{param[:-1]}_type.member_classes",
                f"{param[:-1]}_type.parameters",
            ]
            for param in connectionsList
        }
        columnsDict["sites"] = [
            "name",
            "atom_type.name",
            "atom_type.parameters",
            "charge",
            "mass",
        ]
        if isinstance(columns, list):
            columnsDict = {
                key: columnsDict[key] + columns for key in columnsDict.keys()
            }  # add in any provided columns
    elif format == "specific_columns":
        assert parameters != "all", (
            f"When formatting for specific columns, please set parameter argument to be one of {['sites']+connectionsList}."
            "Otherwise use a format of default."
        )
        columnsDict = {parameter: columns for parameter in parameters}
    elif format == "publication":
        columnsDict = {
            param: [
                "name",
                f"{param[:-1]}_type.member_classes",
                f"{param[:-1]}_type.parameters",
            ]
            for param in connectionsList
        }
        columnsDict["sites"] = [
            "name",
            "atom_type.name",
            "atom_type.parameters",
            "charge",
            "mass",
        ]
        remove_duplicatesBool = True
    elif format == "remove_duplicates":
        assert parameters != "all", (
            f"When formatting for specific columns, please set parameter argument to be one of {['sites']+connectionsList}."
            "Otherwise use a format of default."
        )
        if not columns and parameters == "sites":  # default values
            columns = ["atom_type.name"]
        elif not columns and parameters in connectionsList:
            columns = [f"{parameters[:-1]}_type.member_classes"]
        columnsDict = {parameters: columns}
    else:
        raise ValueError(
            f"Available options for format are 'default', 'specific_columns', 'publication', or 'remove_duplicates'. The incorrect argument passed was {format=}."
        )

    if parameters == "all":
        parametersList = ["sites"] + connectionsList
    elif parameters in connectionsList or parameters == "sites":
        parametersList = [parameters]
    elif isinstance(parameters, list) and all(
        [parameter in connectionsList for parameter in parameters]
    ):
        parametersList = parameters
    else:
        raise ValueError(
            f"parameters argument {parameters} must be one of: 'all', 'sites', {', '.join(connectionsList)}."
        )

    for param in parametersList:
        if not getattr(topology, f"n_{param}"):
            continue
        dataList, columns = _generate_component_lists(
            topology, param, columnsDict.get(param)
        )
        # handle unyts in values
        dataList, columns = _parse_unyts(convert_unyts, dataList, columns)
        dataDict = {col: data for col, data in zip(columns, dataList)}
        outDict[param] = pd.DataFrame(dataDict)  # create dataframe

    if (
        remove_duplicatesBool and topology.n_sites > 0
    ):  # use flag to remove duplicates in sites
        outDict["sites"] = _add_duplicate_indices_to_sites_dataframe(
            outDict["sites"]
        )
        for param in connectionsList:
            if not getattr(topology, f"n_{param}"):
                continue
            outDict[param] = _remove_duplicate_connections(
                outDict[param], param
            )

    if format == "remove_duplicates":
        for df in outDict.values():  # remove duplicate values
            df.drop("Atom Indices", errors="ignore")
            df.drop_duplicates(inplace=True, ignore_index=True)

    return outDict


def _parse_unyts(convert_unyts, dataList, columnsList):
    if convert_unyts == "to_headers":  # move units to the header
        columnsList = _parse_unyts_to_headers(dataList, columnsList)
        dataList = _parse_unyts_no_unytss(dataList)
    elif convert_unyts == "with_data":  # leave units where they are
        pass
    elif convert_unyts == "no_unyts":  # convert units to floats
        dataList = _parse_unyts_no_unytss(dataList)
    else:
        raise ValueError(
            f"Supplied the argument {convert_unyts=} of {type(convert_unyts)}, but must provide one of the arguments 'to_headers', 'with_data', or 'no_unyts'."
        )
    return dataList, columnsList


def _parse_unyts_no_unytss(dataList) -> list:
    for i in range(len(dataList)):
        if isinstance(dataList[i][0], u.unyt_array):
            dataList[i] = [float(x) for x in dataList[i]]  # turn to float
    return dataList


def _parse_unyts_to_headers(dataList, columns) -> list:
    new_colsList = []
    for data, col in zip(dataList, columns):
        if isinstance(data[0], u.unyt_array):
            unit = str(
                data[0].units
            )  # assumption that all data in List is same units
            new_colsList.append(col + f" ({unit})")
        else:
            new_colsList.append(col)
    return new_colsList


def _generate_component_lists(topology, parameter, columns) -> list:
    outList = []
    columnsList = []
    for column in columns:
        valuesList = _recursive_getattr(topology, parameter, column)
        if isinstance(valuesList[0], dict):
            # add keys to columnsList and values to outList
            keys = list(valuesList[0].keys())
            values_dictList = [
                [value[key] for value in valuesList] for key in keys
            ]
            outList.extend(values_dictList)
            columnsList.extend(keys)
        elif isinstance(valuesList[0], u.unyt_array) and not isinstance(
            valuesList[0], u.unyt_quantity
        ):
            outList.extend(np.array(valuesList).T)
            if column == "position":
                columnsList.extend(["x", "y", "z"])
            else:
                columnsList.extend(
                    [f"{column}-({i})" for i in range(len(valuesList[0]))]
                )
        elif isinstance(valuesList[0], tuple) or isinstance(
            valuesList[0], list
        ):  # could be connection_members
            outList.extend(np.array(valuesList).T)
            if "connection_members" in column:
                columnsList.extend(
                    [
                        f"{parameter} member ({i})"
                        for i in range(len(valuesList[0]))
                    ]
                )
            else:
                columnsList.extend(
                    [f"{column}-({i})" for i in range(len(valuesList[0]))]
                )

            # handle positions?
            # handle connection_members
            pass
        else:
            outList.append(valuesList)
            columnsList.append(column)
    return outList, columnsList


def _recursive_getattr(topology, attr, attr_attr):
    """Parse a topology to get a list of attributes from an iterable."""

    def _getattr(obj, attr1):
        try:
            return getattr(obj, attr1)
        except AttributeError as e:
            raise AttributeError(
                f"The GMSO Topology is missing the requested attribute {attr1} from {obj}.{attr_attr}"
            )

    iteritems = getattr(topology, attr)
    parseFunction = lambda x: functools.reduce(
        _getattr, [x] + attr_attr.split(".")
    )
    return list(map(parseFunction, iteritems))


def _add_duplicate_indices_to_sites_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    unique_col = "atom_type.name"  # use to grab what is considered `unique`, may be able to make this a variable in the future
    df["Atom Indices"] = df[unique_col].apply(
        lambda x: ", ".join(
            str(v) for v in df.index[df[unique_col] == x].to_list()
        )
    )
    keep = df[~df.duplicated(subset=unique_col)]
    return keep.reset_index()


def _remove_duplicate_connections(df: pd.DataFrame, parameter) -> pd.DataFrame:
    # dset connection members length
    membersMap = {"bonds": 2, "angles": 3, "dihedrals": 4, "impropers": 4}
    # drop duplicate rows in df
    n_atoms = membersMap[parameter]
    df = df.drop_duplicates(
        subset=[
            f"{parameter[:-1]}_type.member_classes-({i})"
            for i in range(n_atoms)
        ]
    )
    # remove columns for indexing
    # df = df.drop(labels=[f"Atom{i}" for i in range(n_atoms)], axis=1)
    return df.reset_index(drop=True)


def multi_topology_dataframe(topologies: list) -> pd.DataFrame:
    """Take an iterable of topologies and create a combined dataframe to encompass all parameters."""
    assert isinstance(topologies, Iterable)
    assert isinstance(next(iter(topologies)), Topology)
    topList = list(topologies)
    dictList = []
    for top in topList:
        dictList.append(to_dataframeDict(top, format="publication"))
    concatDict = {}
    for parameter in ["sites", "bonds", "angles", "dihedrals", "impropers"]:
        dfsList = list(map(lambda x: x.get(parameter), dictList))
        if not any(elem is not None for elem in dfsList):
            continue
        dfout = pd.concat(
            [
                df.drop("Atom Indices", errors="ignore")
                for df in dfsList
                if not df is None
            ]
        )  # remove missing dfs
        # remove duplicates
        concatDict[parameter] = dfout.drop_duplicates().reset_index()

    return concatDict


def generate_topology_report(topologies: list | Topology) -> pd.DataFrame:
    """Generate information of 2D structure and parameters for an iterable of Topologies."""
