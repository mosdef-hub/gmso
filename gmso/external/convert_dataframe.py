"""Module support for converting to/from Pandas DataFrame objects."""

import copy
import warnings
from operator import attrgetter, itemgetter
import functools
from collections.abc import Iterable

import numpy as np
import unyt as u
from symengine import expand
import pandas as pd

from gmso import Topology
from gmso.core.element import element_by_atomic_number, element_by_symbol
from gmso.core.views import PotentialFilters, get_parameters
from gmso.abc.serialization_utils import unyt_to_dict

pfilter = PotentialFilters.UNIQUE_PARAMETERS
from gmso.exceptions import GMSOError
from gmso.lib.potential_templates import PotentialTemplateLibrary
from gmso.utils.io import import_

pd = import_("pandas")

def to_dataframeDict(topology: Topology, parameter: str="all", format: str="default", columns: list[str]=None, handle_unyts: str="in_headers") -> pd.DataFrame:
    """Return a dictioanry of pandas dataframe objects for a topology.

    Parameters
    ----------
    topology : gmso.Topology, required
        Topology to use for converting values
    parameter : str, optional, default='all'
        A string determining what aspects of the gmso topology will be reported.
        Options are: 'all', 'sites', 'bonds', 'angles', 'dihedrals', and 'impropers'. Defaults to 'all'.
    format : str, optional, default='default'
        The output formatting style for the dataframe. 
        Options are 'default', 'specific_columns', 'publication', `remove_duplicates`. Defaults to 'default'
        'default' will output default column values of ["name", "atom_type.name", "atom_type.parameters", "charge", "mass"],
        and any additional attributes in the `columns` argument.
        'specific_columns' will only output the attributes from the `columns` argument.
        'publication' will use the default outputs, but remove duplicate values from the dataframes. It adds a column labeled
        'Atom Indices' to the `sites` dataframe, which enumerates the indices that the atom_type is a part of.
        `remove_duplicates` will use the labels in passed through the columns argument, and remove duplicates rows in the dataframe.
        You must also set the `parameter` argument to be one of {"sites", "bonds", "angles", "dihedrals", "impropers"}, not {"all"}
        See Notes for more details on what this looks like.
    columns : list of str, optional, default=None
            List of strings that are attributes of the topology site and can be included as entries in the pandas dataframe.
        Examples of these can be found by printing `topology.sites[0].__dict__` or `topology.bonds[0].__dict__`.
        See https://gmso.mosdef.org/en/stable/data_structures.html#gmso.Atom for additional information on labeling.
    handle_unyts: str, optional, default='in_headers'
        The placement/recording of unyt quantities in dataframe. 
        Options are 'in_headers', 'with_data', 'all_floats'
        Determines if numerical values in the DataFrame are saved as unyt quantities or floats. Default method 'in_headers"
        is convert the values to floats, and place a string of the units in the column headers.
        See https://unyt.readthedocs.io/en/stable/usage.html
        for more information about manipulating unyt quantities.

    Returns
    -------
    Dictionary of Pandas Dataframe
        A python dictionary of pandas.Dataframe object, see https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html
        for further information. The keys of this dictionary are the attributes of the topology that are associated with each DataFrame. These
        can be `sites`, `bonds`, `angles`, `dihedrals`, `impropers`, which are determined from the argument `parameter`.

        
    # TODO: Show more examples
    Examples
    ________
    >>> topology.to_dataframe(parameter = 'sites', site_attrs = ['charge'])
        This will return a dataframe with a listing of the sites and include the charges that correspond to each site.
    >>> topology.to_dataframe(parameter = 'dihedrals', site_attrs = ['positions'])
        This will return a dataframe with a listing of the sites that make up each dihedral, the positions of each of
        those sites, and the parameters that are associated with the dihedrals.

    Notes
    ____
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
    outDict = {} # dictionary of dataframes to write out

    # Get columns from format methods
    columnsDict = {}
    connectionsList = ["bonds", "angles", "dihedrals", "impropers"] # these can be handled generally
    remove_duplicatesBool = False # flag to remove duplicate parameters from the dataframe and put indices into a new column
    if format == "default":
        columnsDict = {param: ["name", f"{param[:-1]}_type.member_classes", f"{param[:-1]}_type.parameters"] for param in connectionsList}
        columnsDict["sites"] = ["name", "atom_type.name", "atom_type.parameters", "charge", "mass"]
        if isinstance(columns, list):
            columnsDict = {key:columnsDict[key]+columns for key in columnsDict.keys()} # add in any provided columns
    elif format == "specific_columns":
        assert parameter != "all", (
            f"When formatting for specific columns, please set parameter argument to be one of {['sites']+connectionsList}."
            "Otherwise use a format of default."
        )
        columnsDict = {parameter:columns}
    elif format == "publication":
        columnsDict = {param: ["name", f"{param[:-1]}_type.member_classes", f"{param[:-1]}_type.parameters"] for param in connectionsList}
        columnsDict["sites"] = ["name", "atom_type.name", "atom_type.parameters", "charge", "mass"]
        remove_duplicatesBool = True
    elif format == "remove_duplicates":
        assert parameter != "all", (
            f"When formatting for specific columns, please set parameter argument to be one of {['sites']+connectionsList}."
            "Otherwise use a format of default."
        )
        if not columns and parameter == "sites": # default values
            columns = ["atom_type.name"]
        elif not columns and parameter in connectionsList:
            columns = [f"{parameter[:-1]}_type.member_classes"]
        columnsDict = {parameter:columns}
    else: 
        raise ValueError(f"Please provide formt=['default', 'specific_columns', 'publication']")

    if parameter == "all":
        parametersList = ["sites"] + connectionsList
    elif parameter in connectionsList or parameter == "sites":
        parametersList = [parameter]
    else:
        raise ValueError(f"parameter {parameter} must be one of: 'all', 'sites', {', '.join(connectionsList)}.")

    for param in parametersList:
        if not getattr(topology, f"n_{param}"):
            continue
        dataList, columns = _generate_component_lists(topology, param, columnsDict.get(param))
        # handle unyts in values
        dataList, columns = _parse_unyts(handle_unyts, dataList, columns)
        dataDict = {col:data for col, data in zip(columns, dataList)}
        outDict[param] = pd.DataFrame(dataDict) # create dataframe

    if remove_duplicatesBool and topology.n_sites>0: # use flag to remove duplicates in sites
        outDict["sites"] = _add_duplicate_indices_to_sites_dataframe(outDict["sites"])
        for param in connectionsList:
            if not getattr(topology, f"n_{param}"):
                continue
            outDict[param] = _remove_duplicate_connections(outDict[param], param)
    
    if format == "remove_duplicates":
        for df in outDict.values(): # remove duplicate values
            df.drop("Atom Indices", errors="ignore")
            df.drop_duplicates(inplace=True, ignore_index=True)

    ###############
    # END OF FUNCTION
    ###############
    return outDict
        
def _parse_unyts(handle_unyts, dataList, columnsList):    
    if handle_unyts == "in_headers": # move units to the header
        columnsList = _parse_unyts_to_headers(dataList, columnsList)
        dataList = _parse_unyts_to_floats(dataList)
    elif handle_unyts == "with_data": # leave units where they are
        pass
    elif handle_unyts == "all_floats": # convert units to floats
        dataList = _parse_unyts_to_floats(dataList)
    else:
        raise ValueError(f"Supplied the argument handle_unyts={handle_unyts} of {type(handle_unyts)}, but must provide one of: 'in_headers', 'with_data', or 'all_floats'.")
    return dataList, columnsList

def _parse_unyts_to_floats(dataList) -> list:
    for i in range(len(dataList)):
        if isinstance(dataList[i][0], u.unyt_array):
            dataList[i] = [float(x) for x in dataList[i]] # turn to float
    return dataList

def _parse_unyts_to_headers(dataList, columns) -> list:
    new_colsList = []
    for data, col in zip(dataList, columns):
        if isinstance(data, u.unyt_array):
            unit = str(data[0].units) # assumption that all data in List is same units
            new_colsList.append(col+f"({unit})")
        else:
            new_colsList.append(col)
    return new_colsList

def _generate_component_lists(topology, parameter, columns)-> list:
    outList = []
    columnsList = []
    for column in columns:
        valuesList = _recursive_getattr(topology, parameter, column)
        if isinstance(valuesList[0], dict):
            # add keys to columnsList and values to outList
            keys = list(valuesList[0].keys())
            values_dictList = [[value[key] for value in valuesList] for key in keys]
            outList.extend(values_dictList)
            columnsList.extend(keys)
        elif (isinstance(valuesList[0], u.unyt_array) and not isinstance(valuesList[0], u.unyt_quantity)):
            outList.extend(np.array(valuesList).T)
            if column == "position":
                columnsList.extend(["x", "y", "z"])
            else:
                columnsList.extend([f"{column}-({i})" for i in range(len(valuesList[0]))])
        elif isinstance(valuesList[0], tuple) or isinstance(valuesList[0], list): # could be connection_members
            outList.extend(np.array(valuesList).T)
            if "connection_members" in column:
                columnsList.extend([f"{parameter} member ({i})" for i in range(len(valuesList[0]))])
            else:
                columnsList.extend([f"{column}-({i})" for i in range(len(valuesList[0]))])

            # handle positions?
            # handle connection_members
            pass
        else:
            outList.append(valuesList)
            columnsList.append(column)
    return outList, columnsList

def _recursive_getattr(topology, attr, attr_attr):
    """Parse a topology to get a list of attributes from an iterable"""
    def _getattr(obj, attr1):
        try:
            return getattr(obj, attr1)
        except AttributeError as e:
            raise AttributeError(f"The GMSO Topology is missing the requested attribute {attr1} from {obj}.{attr_attr}")

    iteritems = getattr(topology, attr)
    parseFunction = lambda x: functools.reduce(_getattr, [x] + attr_attr.split('.'))
    return list(map(parseFunction, iteritems))

def _pandas_from_parameters(
    topology, df, parameter, site_attrs=None, unyts_bool=True
    ):
    """Add to a pandas dataframe the site indices for each connection member in a
    multimember topology attribute such as a bond. Also include information about
    those sites in the site_attrs list"""
    if site_attrs is None:
        site_attrs = []
    sites_per_connection = len(
        getattr(self, parameter)[0].connection_members
    )
    for site_index in np.arange(sites_per_connection):
        df["Atom" + str(site_index)] = list(
            str(connection.connection_members[site_index].name)
            + f"({self.get_index(connection.connection_members[site_index])})"
            for connection in getattr(self, parameter)
        )
    for attr in site_attrs:
        df = self._parse_dataframe_attrs(
            df, attr, parameter, sites_per_connection, unyts_bool
        )
    return df

def _parse_dataframe_attrs(
    self, df, attr, parameter, sites_per_connection=1, unyts_bool=True
    ):
    """Parses an attribute string to correctly format and return the topology attribute
    into a pandas dataframe"""
    if parameter == "sites":
        if "." in attr:
            try:
                attr1, attr2 = attr.split(".")
                df[attr] = list(
                    _return_float_for_unyt(
                        getattr(getattr(site, attr1), attr2),
                        unyts_bool,
                    )
                    for site in self.sites
                )
            except AttributeError:
                raise AttributeError(
                    f"The attribute {attr} is not in this gmso object."
                )
        elif attr == "positions" or attr == "position":
            for i, dimension in enumerate(["x", "y", "z"]):
                df[dimension] = list(
                    _return_float_for_unyt(
                        getattr(site, "position")[i], unyts_bool
                    )
                    for site in self.sites
                )
        elif attr == "charge" or attr == "charges":
            df["charge (e)"] = list(
                site.charge.in_units(
                    u.Unit(
                        "elementary_charge", registry=UnitReg.default_reg()
                    )
                ).to_value()
                for site in self.sites
            )
        else:
            try:
                df[attr] = list(
                    _return_float_for_unyt(getattr(site, attr), unyts_bool)
                    for site in self.sites
                )
            except AttributeError:
                raise AttributeError(
                    f"The attribute {attr} is not in this gmso object."
                )

    elif parameter in ["bonds", "angles", "dihedrals", "impropers"]:
        for site_index in np.arange(sites_per_connection):
            if "." in attr:
                try:
                    attr1, attr2 = attr.split(".")
                    df[attr + " Atom" + str(site_index)] = list(
                        _return_float_for_unyt(
                            getattr(
                                getattr(
                                    connection.connection_members[
                                        site_index
                                    ],
                                    attr1,
                                ),
                                attr2,
                            ),
                            unyts_bool,
                        )
                        for connection in getattr(self, parameter)
                    )
                except AttributeError:
                    raise AttributeError(
                        f"The attribute {attr} is not in this gmso object."
                    )
            elif attr == "positions" or attr == "position":
                df["x Atom" + str(site_index) + " (nm)"] = list(
                    _return_float_for_unyt(
                        getattr(
                            connection.connection_members[site_index],
                            "position",
                        )[0],
                        unyts_bool,
                    )
                    for connection in getattr(self, parameter)
                )
                df["y Atom" + str(site_index) + " (nm)"] = list(
                    _return_float_for_unyt(
                        getattr(
                            connection.connection_members[site_index],
                            "position",
                        )[1],
                        unyts_bool,
                    )
                    for connection in getattr(self, parameter)
                )
                df["z Atom" + str(site_index) + " (nm)"] = list(
                    _return_float_for_unyt(
                        getattr(
                            connection.connection_members[site_index],
                            "position",
                        )[2],
                        unyts_bool,
                    )
                    for connection in getattr(self, parameter)
                )
            elif attr == "charge" or attr == "charges":
                df["charge Atom" + str(site_index) + " (e)"] = list(
                    getattr(
                        connection.connection_members[site_index],
                        "charge",
                    )
                    .in_units(
                        u.Unit(
                            "elementary_charge",
                            registry=UnitReg.default_reg(),
                        )
                    )
                    .value
                    for connection in getattr(self, parameter)
                )
            else:
                try:
                    df[f"{attr} Atom {site_index}"] = list(
                        _return_float_for_unyt(
                            getattr(
                                connection.connection_members[site_index],
                                attr,
                            ),
                            unyts_bool,
                        )
                        for connection in getattr(self, parameter)
                    )
                except AttributeError:
                    raise AttributeError(
                        f"The attribute {attr} is not in this gmso object."
                    )
    else:
        raise AttributeError(
            f"{parameter} is not yet supported for adding labels to a dataframe. \
                Please use  one of 'sites', 'bonds', 'angles', 'dihedrals', or 'impropers'"
        )
    return df

def _parse_parameter_expression(self, df, parameter, unyts_bool):
    """Take a given topology attribute and return the parameters associated with it"""
    for i, param in enumerate(
        getattr(
            getattr(self, parameter)[0], parameter[:-1] + "_type"
        ).parameters
    ):
        df[
            f"Parameter {i} ({param}) {getattr(getattr(self, parameter)[0], parameter[:-1]+'_type').parameters[param].units}"
        ] = list(
            _return_float_for_unyt(
                getattr(connection, parameter[:-1] + "_type").parameters[
                    param
                ],
                unyts_bool,
            )
            for connection in getattr(self, parameter)
        )
    return df

def _return_float_for_unyt(unyt_quant, unyts_bool):
    try:
        return unyt_quant if unyts_bool else unyt_to_dict(unyt_quant)["array"]
    except TypeError:
        return unyt_quant
    
def _add_duplicate_indices_to_sites_dataframe(df:pd.DataFrame) -> pd.DataFrame:
    unique_col = "atom_type.name" # use to grab what is considered `unique`, may be able to make this a variable in the future
    df["Atom Indices"] = df[unique_col].apply(lambda x: ", ".join(str(v) for v in df.index[df[unique_col] == x].to_list()))
    keep = df[~df.duplicated(subset=unique_col)]
    return keep.reset_index()

def _remove_duplicate_connections(df:pd.DataFrame, parameter) -> pd.DataFrame:
    # dset connection members length
    membersMap = {
        "bonds":2, "angles":3, "dihedrals":4, "impropers":4
    }
    # drop duplicate rows in df
    n_atoms = membersMap[parameter]
    df = df.drop_duplicates(subset=[f"{parameter[:-1]}_type.member_classes-({i})" for i in range(n_atoms)])
    # remove columns for indexing
    #df = df.drop(labels=[f"Atom{i}" for i in range(n_atoms)], axis=1)
    return df.reset_index(drop=True)
    
def multi_topology_dataframe(topologies : list) -> pd.DataFrame:
    """Take an iterable of topologies and create a combined dataframe to encompass all parameters."""
    assert isinstance(topologies, Iterable)
    assert isinstance(next(iter(topologies)), Topology)
    topList = list(topologies)
    dictList = []
    for top in topList:
        dictList.append(to_dataframeDict(top, format="publication"))
    concatDict = {}
    for parameter in ["sites","bonds", "angles", "dihedrals", "impropers"]:
        dfsList = list(map(lambda x: x.get(parameter), dictList))
        if not any(elem is not None for elem in dfsList):
            continue
        dfout = pd.concat([df.drop("Atom Indices", errors="ignore") for df in dfsList if not df is None]) # remove missing dfs
        # remove duplicates
        concatDict[parameter] = dfout.drop_duplicates().reset_index()

    return concatDict


def generate_topology_report(topologies : list | Topology) -> pd.DataFrame:
    """Generate information of 2D structure and parameters for an iterable of Topologies."""
    pass
    