================================
DiffDock scipion plugin
================================

**Documentation under development, sorry for the inconvenience**

Scipion framework plugin for the use of tools provided by DiffDock software

===================
Install this plugin
===================

You will need to use `3.0.0 <https://github.com/I2PC/scipion/releases/tag/v3.0>`_ version of Scipion
to run these protocols. To install the plugin, you have two options:

- **Stable version**  

.. code-block:: 

      scipion installp -p scipion-chem-diffdock
      
OR

  - through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**
      
- **Developer's version** 

1. **Download repository**:

.. code-block::

            git clone https://github.com/scipion-chem/scipion-chem-diffdock.git

2. **Switch to the desired branch** (main or devel):

Scipion-chem-diffdock is constantly under development.
If you want a relatively older an more stable version, use main branch (default).
If you want the latest changes and developments, user devel branch.

.. code-block::

            cd scipion-chem-diffdock
            git checkout devel

3. **Install**:

.. code-block::

            scipion installp -p path_to_scipion-chem-diffdock --devel

- **Binary files** 

Atom_struct_utils plugin is a pure Python module, no binary files are required. 

- **Tests**

To check the installation, simply run the following Scipion test:

===============
Buildbot status
===============

Status devel version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/bioinformatics_dev.svg

Status production version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/bioinformatics_prod.svg
