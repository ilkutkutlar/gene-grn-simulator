from PyQt5.QtWidgets import QWidget, QFormLayout, QLineEdit, QLabel, QDialog, QPushButton

"""
A simple input dialog which allows multiple 
"""


class QMultipleInputDialog(QDialog):

    def __init__(self, labels, handler, button_text="OK"):
        super().__init__()

        self.handler = handler
        (self.main, self.edits) = self._make_form(labels)
        self.main.addWidget(self._make_button(button_text, handler))

        self.setLayout(self.main)
        self.show()

    def _make_form(self, labels):
        form = QFormLayout()
        edits = []

        for l in labels:
            label = QLabel(l)
            edit = QLineEdit()
            edits.append(edit)
            form.addRow(label, edit)

        return form, edits

    def _make_button(self, button_text, click_handler):
        ok_button = QPushButton(button_text)
        ok_button.clicked.connect(self.ok_button_clicked)
        return ok_button

    def ok_button_clicked(self):
        inputs = [str(e.text()) for e in self.edits]
        self.handler(inputs)
        self.close()
