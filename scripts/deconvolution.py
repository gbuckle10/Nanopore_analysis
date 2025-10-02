import logging

class Deconvolution:
    def __init__(self, config, tool_paths):
        """
        Initialise deconvolution handler.
        """
        self.config = config
        self.tool_paths = tool_paths
        self.logger = logging.getLogger('pipeline')