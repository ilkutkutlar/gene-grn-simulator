from PyQt5.QtWidgets import QWidget, QFormLayout, QLineEdit, QLabel, QDialog


class QMultipleInputDialog(QDialog):

    def __init__(self, labels):
        super().__init__()
        self.main_layout = QFormLayout()

        for l in labels:
            label = QLabel()
            label.setText(l)
            edit = QLineEdit()
            self.main_layout.addRow(label, edit)

        self.setLayout(self.main_layout)
        self.show()
