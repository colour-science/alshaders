import pymel.core as pm
from alShaders import alShadersTemplate

class AEalCombineFloatTemplate(alShadersTemplate):
    def setup(self):
        self.addSwatch()
        self.beginScrollLayout()

        self.addControl('input1', label='Input 1')
        self.addControl('input2', label='Input 2')
        self.addControl('input3', label='Input 3')
        self.addControl('combineOp', label='Operation')

        pm.mel.AEdependNodeTemplate(self.nodeName)

        self.addExtraControls()
        self.endScrollLayout()