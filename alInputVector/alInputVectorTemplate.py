import pymel.core as pm
from alShaders import alShadersTemplate

class AEalInputVectorTemplate(alShadersTemplate):
    def setup(self):
        self.addSwatch()
        self.beginScrollLayout()

        self.addControl('input', label='Input')
        self.addControl('userName', label='User name')
        self.addControl('vector', label='Custom')
        self.addControl('type', label='Type')

        self.beginLayout('Transform', collapse=True)
        self.addControl('matrix', label='Matrix')
        self.addControl('coordinates')
        self.endLayout() # end Transform

        pm.mel.AEdependNodeTemplate(self.nodeName)

        self.addExtraControls()
        self.endScrollLayout()