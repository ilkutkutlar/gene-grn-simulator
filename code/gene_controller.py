from models import Network


class GeneController:
    instance = None

    def __init__(self):
        self.network = Network()

    @staticmethod
    def get_instance():
        if not GeneController.instance:
            GeneController.instance = GeneController()
        return GeneController.instance
