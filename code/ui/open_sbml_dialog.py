from PyQt5.QtWidgets import QFileDialog, QMessageBox

from input_output.sbml_parser import SbmlParser
from ui.gene_controller import GeneController


class OpenSbmlDialog(QFileDialog):

    def __init__(self, parent):
        super().__init__()

        self.setFileMode(QFileDialog.AnyFile)

        # TODO: Why is the filter not working?
        # file_dialog.setFilter("SBML files (*.xml)")

        if self.exec_():
            # selectedFiles returns all files, we only want to open one
            filename: str = self.selectedFiles()[0]
            net = SbmlParser.parse(filename)
            GeneController.get_instance().network = net

            # Display a nice message showing the contents of the file loaded
            message = QMessageBox()
            message.setText("SBML file has been opened and this network has been loaded: \n\n" + net.__str__())
            message.exec_()

            parent.species_tab.update_list()
            parent.reactions_tab.update_list()
