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
