.. PCMDI Metrics Package documentation master file, created by
   sphinx-quickstart on Thu May 24 11:23:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

The PCMDI Metrics Package (PMP)
=================================================

Overview - science
------------------

The PCMDI metrics package (PMP) is designed to objectively compare results from Earth System Models (ESMs) with observations.    
It includes an increasing diverse set of relatively robust high level summary statistics that gauge a wide range of model behaviour.   
The focus on the PMP is on the physical characteristics of ESM behavior.

PCMDI scientists use the PMP to synthesize model behaivor (e.g.,doi:10.1007/s00382-018-4355-4). PCMDI also collaborates with community-based expert teams to help identify what analysis to include in the PMP.  The PCMDI uses the PMP to produce and make summaries publically available (https://pcmdi.llnl.gov/research/metrics).  PCMDI provides support to modeling groups who use the PMP in their model development process.  

PCMDI's simulation summaries are produced in the context of all model simulations contributed to CMIP6 and earlier CMIP phases. Among other purposes, this enables modeling groups to evaluate changes during the development cycle in the context of the structural error distribution of the multi-model ensemble. Currently, the comparisons emphasize large- to global-scale annual cycle performance metrics. Current work in v1.x development branches include established statistics for ENSO, regional monsoon precipitation, and the diurnal cycle of precipitation. These diagnostics will be included in a future PMP release.

The metrics package consists of four parts: 1) Analysis software, 2) an observationally-based database of global (or near global, land or ocean) annual cycle climatologies, 3) a database of performance metrics computed for CMIP models and 4) package documentation. The package expects model data to be CF-compliant. To successfully use the package some input data "conditioning" may be required. We provide several demo scripts within the package to facilitate new users.

Overview - software
-------------------

The PMP consists of four parts: 1) Analysis software, 2) an observationally-based database of global (or near global, land or ocean) annual cycle climatologies, 3) a database of performance metrics computed for CMIP models and 4) package documentation. The package expects model data to be CF-compliant. To successfully use the package some input data "conditioning" may be required. We provide several demo scripts to facilitate new users.




Current State
-------------
Users of the current release (v1.1.2) will still need to contact the PMP developers (pcmdi-metrics@llnl.gov) to obtain supporting datasets and get started using the package.

v1.1.2 - Includes high frequency precipitation characteristics including the diurnal cycle and 'intermittency'.

Earlier versions include:

v1.1.2 - Now managed through Anaconda, and tied to UV-CDAT 2.10. Weights on bias statistic added. Extensive provenance information incorporated into json files.

v1.1 - First public release, emphasizing climatological statistics, with development branches for ENSO and regional monsoon precipitation indices

v1.0 - Prototype version of the PMP


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   example.md


API
===
.. toctree::
   :maxdepth: 0

   API/api.rst

Jupyter Notebooks
=================
.. toctree::
  :maxdepth: 0

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
