import pymel.core as pm
from alShaders import alShadersTemplate

class AEalHairTemplate(alShadersTemplate):
    def setup(self):
        self.addSwatch()
        self.beginScrollLayout()

        #self.beginLayout('Diffuse', collapse=False)
        #self.addControl('diffuseStrength', label='Strengh')
        #self.addControl('diffuseColor', label='Color')
        #self.endLayout() # end Diffuse

        self.beginLayout('Fibre properties', collapse=False)
        self.addControl('hairColor', label='Color')
        self.addControl('hairColorDensity', label='Color density')
        self.addControl('specularWidth', label='Highlight width')
        self.addControl('specularShift', label='Highlight shift')
        self.endLayout() # end Fibre properties

        self.beginLayout('Specular 1', collapse=False)
        self.addControl('specular1Strength', label='Strength')
        self.addControl('specular1Color', label='Tint')
        self.endLayout() # end Specular 1

        self.beginLayout('Specular 2', collapse=False)
        self.addControl('specular2Strength', label='Strength')
        self.addControl('specular2Color', label='Tint')
        self.addControl('glintStrength', label='Glint strength')
        self.addControl('glintRolloff', label='Glint rolloff')
        self.endLayout() # end Specular 2

        self.beginLayout('Transmission', collapse=False)
        self.addControl('transmissionStrength', label='Strength')
        self.addControl('transmissionColor', label='Tint')
        self.addControl('transmissionRolloff', label='Rolloff')
        self.endLayout() # end Transmission

        self.beginLayout('Dual scattering', collapse=False)
        self.addControl('densityFront')
        self.addControl('densityBack')
        self.endLayout() # end Dual scattering

        self.beginLayout('Advanced', collapse=True)
        self.addControl('extraSamples', label='Extra samples')
        self.endLayout() # end Advacned

        pm.mel.AEdependNodeTemplate(self.nodeName)

        self.addExtraControls()
        self.endScrollLayout()