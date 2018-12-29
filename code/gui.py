import sys

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget, QApplication, QLabel, QVBoxLayout, QListWidget, QHBoxLayout, QPushButton, \
    QInputDialog, QDialog, QComboBox, QFormLayout, QLineEdit, QMainWindow

from gene_controller import GeneController
from models import TranslationReaction, TranscriptionReaction, MrnaDegradationReaction, ProteinDegradationReaction, \
    CustomReaction


class AddReactionDialog(QDialog):
    def __init__(self):
        super().__init__()
        self._init_ui()

    def _init_combo(self):
        self.combo.addItem("Transcription Reaction")
        self.combo.addItem("Translation Reaction")
        self.combo.addItem("mRNA Degradation Reaction")
        self.combo.addItem("Protein Degradation Reaction")
        self.combo.addItem("Custom Reaction")
        self.combo.currentIndexChanged.connect(self._handler_reaction_type_changed)

    def _init_form_fields(self):
        self.fields = {}

        # TODO: Regulation!
        # Common
        self.fields["name"] = QLineEdit()
        self.fields["name"].setPlaceholderText("Name")

        self.fields["rp_info"] = QLabel()
        self.fields["rp_info"].setText("Reactants and products must be comma separated names of species")

        self.fields["reactants"] = QLineEdit()
        self.fields["reactants"].setPlaceholderText("Reactants")

        self.fields["products"] = QLineEdit()
        self.fields["products"].setPlaceholderText("Products")

        # Transcription
        self.fields["transcription_rate"] = QLineEdit()
        self.fields["transcription_rate"].setPlaceholderText("Transcription Rate")

        self.fields["kd"] = QLineEdit()
        self.fields["kd"].setPlaceholderText("Kd")

        self.fields["hill_coefficient"] = QLineEdit()
        self.fields["hill_coefficient"].setPlaceholderText("Hill coefficient")

        # Translation
        self.fields["translation_rate"] = QLineEdit()
        self.fields["translation_rate"].setPlaceholderText("Tranlation Rate")

        # mRNA decay
        self.fields["mrna_decay_rate"] = QLineEdit()
        self.fields["mrna_decay_rate"].setPlaceholderText("Decay Rate")

        # Protein decay
        self.fields["protein_decay_rate"] = QLineEdit()
        self.fields["protein_decay_rate"].setPlaceholderText("Decay Rate")

        # Custom reaction
        self.fields["custom_equation"] = QLineEdit()
        self.fields["custom_equation"].setPlaceholderText("Equation")

        # rp_layout = QHBoxLayout()
        # rp_layout.addWidget(self.fields["reactants"])
        # arrow_label = QLabel()
        # arrow_label.setText("->")
        # rp_layout.addWidget(arrow_label)
        # rp_layout.addWidget(self.fields["products"])

        for field in self.fields:
            self.fields[field].setVisible(False)
            self.form.addWidget(self.fields[field])


        # Make the initial fields visible (i.e. the combo box has this reaction
        # selected when the dialog is first opened)
        self.fields["name"].setVisible(True)
        self.fields["reactants"].setVisible(True)
        self.fields["products"].setVisible(True)
        self.fields["transcription_rate"].setVisible(True)
        self.fields["kd"].setVisible(True)
        self.fields["hill_coefficient"].setVisible(True)

    def _handler_reaction_type_changed(self, index):
        for field in self.fields:
            self.fields[field].setVisible(False)

        self.fields["name"].setVisible(True)
        self.fields["reactants"].setVisible(True)
        self.fields["products"].setVisible(True)

        if index == 0:
            self.fields["transcription_rate"].setVisible(True)
            self.fields["kd"].setVisible(True)
            self.fields["hill_coefficient"].setVisible(True)
        elif index == 1:
            self.fields["translation_rate"].setVisible(True)
        elif index == 2:
            self.fields["mrna_decay_rate"].setVisible(True)
        elif index == 3:
            self.fields["protein_decay_rate"].setVisible(True)
        elif index == 4:
            self.fields["custom_equation"].setVisible(True)

    def _handler_ok_button_clicked(self):
        index = self.combo.currentIndex()

        if index == 0:
            GeneController.get_instance().network.reactions.append(
                TranscriptionReaction(self.fields["transcription_rate"].text(),
                                      self.fields["kd"].text(),
                                      self.fields["hill_coefficient"].text(), left=[], right=[])
            )
        elif index == 1:
            GeneController.get_instance().network.reactions.append(
                TranslationReaction(self.fields["translation_rate"].text(), left=[], right=[])
            )
        elif index == 2:
            GeneController.get_instance().network.reactions.append(
                MrnaDegradationReaction(self.fields["mrna_decay_rate"].text(), left=[], right=[])
            )
        elif index == 3:
            GeneController.get_instance().network.reactions.append(
                ProteinDegradationReaction(self.fields["protein_decay_rate"].text(), left=[], right=[])
            )
        elif index == 4:
            GeneController.get_instance().network.reactions.append(
                CustomReaction(self.fields["custom_equation"].text(), left=[], right=[])
            )

        self.close()

    def _init_ui(self):
        self.layout = QVBoxLayout()  # The main layout
        self.form = QVBoxLayout()  # This holds the text fields for reaction parameters
        self.combo = QComboBox()  # Choose the type of reaction
        self.ok_button = QPushButton("OK", self)
        self.ok_button.clicked.connect(self._handler_ok_button_clicked)

        self._init_combo()
        self._init_form_fields()

        self.layout.addWidget(self.combo)
        self.layout.addLayout(self.form)
        self.layout.addWidget(self.ok_button)

        self.setLayout(self.layout)

        self.setWindowTitle("Add reaction")
        self.setWindowModality(Qt.WindowModal)


class GeneWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self._init_ui()

    def _init_ui(self):
        layout = QHBoxLayout()
        layout.addLayout(self._species_panel())
        layout.addLayout(self._reactions_panel())
        self._init_menubar()

        central_widget = QWidget()
        central_widget.setLayout(layout)

        self.setCentralWidget(central_widget)
        self.setWindowTitle("Gene")
        self.show()

    def _init_menubar(self):
        self.menubar = self.menuBar()
        file = self.menubar.addMenu("File")
        file.addAction("Open SBML file")
        simulate = self.menubar.addMenu("Simulate")
        simulate.addAction("Deterministic (ODE)")
        simulate.addAction("Stochastic (Gillespie Algorithm)")

    def _reactions_dialog_returned(self):
        self._refresh_reactions_list()

    def _handler_add_reactions_button(self):
        dialog = AddReactionDialog()
        dialog.finished.connect(self._reactions_dialog_returned)
        dialog.exec_()

    def _handler_add_species_button(self):
        (species, ok) = QInputDialog.getText(self, 'Add new species', 'Species name:')

        if ok:
            GeneController.get_instance().network.species[species] = 0
            self._refresh_species_list()

    def _refresh_reactions_list(self):
        self.reactions_list.clear()
        for reaction in GeneController.get_instance().network.reactions:
            self.reactions_list.addItem(reaction.__str__())

    def _reactions_panel(self):
        self.reactions_label = QLabel(self)
        self.reactions_label.setText("Reactions")
        self.reactions_list = QListWidget()

        self._refresh_reactions_list()

        self.add_button = QPushButton("Add")
        self.remove_button = QPushButton("Remove")
        self.add_button.clicked.connect(self._handler_add_reactions_button)

        self.reactions_layout = QVBoxLayout()
        self.reactions_layout.addWidget(self.reactions_label)
        self.reactions_layout.addWidget(self.reactions_list)
        self.reactions_layout.addWidget(self.add_button)
        self.reactions_layout.addWidget(self.remove_button)
        self.reactions_layout.addStretch()

        return self.reactions_layout

    def _refresh_species_list(self):
        self.species_list.clear()
        for species in GeneController.get_instance().network.species:
            self.species_list.addItem(species)

    def _species_panel(self):
        self.species_label = QLabel(self)
        self.species_label.setText("Species")
        self.species_list = QListWidget()

        self._refresh_species_list()

        self.add_species_button = QPushButton("Add")
        self.remove_species_button = QPushButton("Remove")
        self.add_species_button.clicked.connect(self._handler_add_species_button)

        self.species_layout = QVBoxLayout()
        self.species_layout.addWidget(self.species_label)
        self.species_layout.addWidget(self.species_list)
        self.species_layout.addWidget(self.add_species_button)
        self.species_layout.addWidget(self.remove_species_button)
        self.species_layout.addStretch()
        return self.species_layout


app = QApplication([])
g = GeneWindow()
app.exec_()
