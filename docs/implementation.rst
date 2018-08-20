######################
Implementation details
######################


All functionality is contained in the lib directory:

    * lib/
      Contains general functions and usage elements. It is encouraged to put
      higher complexity functions and special operations not directly related to
      the content here. Files should bundle specific functionalities.

    * lib/api
      Connections with external data sources and functionality. Both file based
      and network based sources should be implemented as separate files. If
      state is required to operate these sources, such as at minimum a filepath,
      it is encouraged to create classes to wrap these sources.

    * lib/model
      Data structures used in the process. Many layered data structures and
      structures not well representing their later usage should be abstracted
      here.

API Components structure
========================

API components are built similarly to a singleton design pattern with global
instances that are per default empty. Instances can be enabled by calling the
configure function. All API classes implementing singleton behavior inherit from
LazyConfigure in :code:`lib/singleton.py`.

This structure enables the shared usage of global API instances, which greatly
reduces the complexity of arguments needed to process case and json level calls
in most common situations, where API configuration is done through either
command-line arguments or configuration file, both of which having global
validity.

In addition to API-components requiring special credentials and paths. Internal
components such as the errorfixer are also implemented as singletons.


jSON Parsing and Case generation
================================

The general structure of the jSON parsing structure and case object abstraction
have been constructed with changing jSON schemas for external information in
mind.

json_parser - jSON-level IO
---------------------------

The *json_parser* module provides classes for the interpretation of raw
json information provided by Face2Gene. All functionality directly related to
quirks in the raw data is to be processed by these classes. Export of case
object information is also handled here.

If new fields are added to the raw json format a new getter should be
implemented or an existing one extended as far the bundling makes sense. Care
should be taken to expose an abstracted interface to the case object. Direct
json level access is not recommended and should be refactored with time.


case - Central model of case information
----------------------------------------

The case object bundles information from a large number of sources to provide
the information needed to generate outputs and information in the preprocessing
steps.

Integration of various data sources and specific shaping of data should be
handled here. For example phenomization combines information from the Face2Gene
detected syndromes with further suggestions from Phenomizer. Both data sources
are input as dataframes and merged on specific keys. Additional transforms are
done to preserve uniform naming and omission of duplicate information.

It is easy to leave too much functionality in the case object, which does not
work well in multiprocessing issues, where serialisability of the case object is
of special concern. Complex functionality should be considered to be implemented
in a general way in a separate module and exposed in the case object via
object-specific methods. For example VCF generation has been separated into a
separate module to reduce code complexity and facilitate testing of individual
components.

hgvs_model - Special functionality for HGVS parsing
---------------------------------------------------

Variant information provided by Face2Gene presents considerable complexity and
has been separated into a separate representation with a lot of field specific
parsing. This component attempts to integrate Mutalyzer information with
manually corrected data and some heuristics for generation of mostly correct
hgvs codes.

Changes to this component should be done with case and should be tested
thoroughly against existing HGVS sequences.

config - Global configuration
-----------------------------

The PEDIAConfig object merges information from the config.ini file with
optionally parsed command line arguments.

New config options should be interpreted here and exposed through an abstracted
interface as instance attributes of the configuration object. Direct access to
the configuration object and commandline arguments is highly discouraged.


Configuration and external usage interface
==========================================

Currently only one usage implementation is availabe in **preprocess.py**. There
exists an partial implementation for the inspection of case information in the
Trello bot implementation. Which parses important information from cases for
visual output in the trello board interface.

Generally an implementation flow involves the following steps.

- Configure needed singleton API objects.
- Create json objects
- Create case objects
- Do operations on case objects
- Export to some other data format
