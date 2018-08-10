from PyQt5 import QtWidgets
from h5xplorer.menu_tools import *
from h5xplorer.menu_plot import *

def context_menu(self, treeview, position):

    """Generate a right-click menu for the items"""

    all_item = get_current_item(self,treeview,single=False)

    if len(all_item) == 1:

        item = all_item[0]
        grp = get_current_hdf5_group(self,item)
        data = get_group_data(grp)
        name = grp.name.strip('/')
        print(name)

        if name == 'PDBsim':
            list_operations = ['Print attrs','-','Connection Graph']

        elif name.startswith('PDBunique/cluster'):
            list_operations = ['Print attrs','-','Protein Graph']

        elif data is not None:

            if data.ndim == 1:
                list_operations = ['Print attrs','-','Plot Hist', 'Plot Line']

            if data.ndim == 2:
                list_operations = ['Print attrs','-','Plot Hist', 'Plot Map']

        else:
            list_operations = ['Print attrs']

        action,actions = get_actions(treeview,position,list_operations)

        if action == actions['Print attrs']:
            send_dict_to_console(self,item,treeview)

        if 'Plot Hist' in actions:
            if action == actions['Plot Hist']:
                plot_histogram(self,item,treeview)

        if 'Plot Line' in actions:
            if action == actions['Plot Line']:
                plot_line(self,item,treeview)

        if 'Plot Map' in actions:
            if action == actions['Plot Map']:
                plot2d(self,item,treeview)

        if 'Connection Graph' in actions:
            if action == actions['Connection Graph']:

                data_dict = {'_grp':grp}
                treeview.emitDict.emit(data_dict)

                cmd = 'connection_graph(_grp)'
                data_dict = {'exec_cmd':cmd}
                treeview.emitDict.emit(data_dict)

        if 'Protein Graph' in actions:
            if action == actions['Protein Graph']:

                data_dict = {'_grp':grp}
                treeview.emitDict.emit(data_dict)

                cmd = 'protein_graph(_grp)'
                data_dict = {'exec_cmd':cmd}
                treeview.emitDict.emit(data_dict)