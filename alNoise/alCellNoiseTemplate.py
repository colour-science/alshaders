import pymel.core as pm
from alShaders import alShadersTemplate

class AEalCellNoiseTemplate(alShadersTemplate):
    def setup(self):
        self.addSwatch()
        self.beginScrollLayout()

        self.addControl('space', label='Space')
        self.addControl('frequency', label='Frequency')
        self.addControl('mode', label='Mode')

        self.beginLayout('Features', collapse=False)
        self.addControl('octaves', label='Octaves')
        self.addControl('lacunarity', label='Lacunarity')
        self.addControl('randomness', label='Randomness')

        self.beginLayout('Feature weights', collapse=False)
        self.addControl('f1w', label='F1')
        self.addControl('f2w', label='F2')
        self.addControl('f3w', label='F3')
        self.addControl('f4w', label='F4')
        self.addControl('mynkowskiShape', label='Shape')
        self.endLayout() # end Feature weights

        self.addRemapControls()

        self.addControl('color1', label='Color 1')
        self.addControl('color2', label='Color 2')

        self.endLayout() # end Features

        self.beginLayout('Chips', collapse=False)
        self.addControl('chipColor1', label='Color 1')
        self.addControl('chipProb1', label='Probability 1')
        self.addControl('chipColor2', label='Color 2')
        self.addControl('chipProb2', label='Probability 2')
        self.addControl('chipColor3', label='Color 3')
        self.addControl('chipProb3', label='Probability 3')
        self.addControl('chipColor4', label='Color 4')
        self.addControl('chipProb4', label='Probability 4')
        self.addControl('chipColor5', label='Color 5')
        self.addControl('chipProb5', label='Probability 5')
        self.endLayout() # end Chips

        self.addControl('P')

        pm.mel.AEdependNodeTemplate(self.nodeName)

        self.addExtraControls()
        self.endScrollLayout()