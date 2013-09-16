import pymel.core as pm
from alShaders import alShadersTemplate

class AEalSurfaceTemplate(alShadersTemplate):
    def setup(self):
        self.addSwatch()
        self.beginScrollLayout()

        self.beginLayout("Matte", collapse=True)
        self.addControl("aiEnableMatte", label="Enable Matte")
        self.addControl("aiMatteColor", label="Matte Color")
        self.addControl("aiMatteColorA", label="Matte Opacity")
        self.endLayout()
        
        self.beginLayout('Diffuse', collapse=False)

        self.addControl('diffuseStrength', label='Strength')
        self.addControl('diffuseColor', label='Color')
        self.addControl('diffuseRoughness', label='Roughness')

        self.beginLayout('Backlight', collapse=True)
        self.addControl('backlightStrength', label='Strength')
        self.addControl('backlightColor', label='Color')
        self.endLayout() # end Backlight
    
        self.beginLayout('Subsurface scattering')
        self.addControl('sssMix', label='Mix')
        self.addControl('sssRadius', label='Distance')
        self.addControl('sssRadiusColor', label='Color')
        self.addControl('sssDensityScale', label='Density scale')
        self.endLayout() # end Subsurface scattering

        self.beginLayout('Advanced')
        self.addControl('diffuseExtraSamples', label='Extra samples')
        self.addControl('diffuseEnableCaustics', label='Enable caustics')
        self.endLayout() # end Advanced

        self.endLayout() # end Diffuse

        self.beginLayout('Specular 1', collapse=False)
        self.addControl('specular1Strength', label='Strength')
        self.addControl('specular1Color', label='Color')
        self.addControl('specular1Roughness', label='Roughness')
        self.addControl('specular1Ior', label='IOR')

        self.beginLayout('Advanced')
        self.addControl('specular1RoughnessDepthScale', label='Roughness depth scale')
        self.addControl('specular1ExtraSamples', label='Extra samples')
        self.addControl('specular1Normal', label='Normal')
        self.endLayout() # end Advanced

        self.endLayout() # end Specular 1

        self.beginLayout('Specular 2')
        self.addControl('specular2Strength', label='Strength')
        self.addControl('specular2Color', label='Color')
        self.addControl('specular2Roughness', label='Roughness')
        self.addControl('specular2Ior', label='IOR')

        self.beginLayout('Advanced')
        self.addControl('specular2RoughnessDepthScale', label='Roughness depth scale')
        self.addControl('specular2ExtraSamples', label='Extra samples')
        self.addControl('specular2Normal', label='Normal')
        self.endLayout() # end Advanced
        
        self.endLayout() # end Specular 1

        self.beginLayout('Transmission')
        self.addControl('transmissionStrength', label='Strength')
        self.addControl('transmissionColor', label='Color')
        self.addControl('transmissionLinkToSpecular1', label='Link to spec')
        self.addControl('transmissionRoughness', label='Roughness')
        self.addControl('transmissionIor', label='IOR')

        self.beginLayout('Scattering')
        self.addControl('ssStrength', label='Strength')
        self.addControl('ssTargetColor', label='Color')
        self.addControl('ssDirection', 'Direction')
        self.addControl('ssBalance', label='Balance')
        self.addControl('ssInScattering', label='In-scattering')
        self.addControl('ssDensityScale', 'Density scale')
        self.addControl('ssSpecifyCoefficients', label='Specify coefficients')
        self.addControl('ssScattering', label='Scattering')
        self.addControl('ssAbsorption', label='Absorption')
        self.endLayout() # end Scattering

        self.beginLayout('Advanced', collapse=True)
        self.addControl('transmissionRoughnessDepthScale', label='Roughness depth scale')
        self.addControl('transmissionExtraSamples', label='Extra samples')
        self.addControl('transmissionEnableCaustics', label='Enable internal reflections')
        self.endLayout() # end Advanced

        self.endLayout() # end Transmission

        self.beginLayout('Emission', collapse=True)
        self.addControl('emissionStrength', label='Strength')
        self.addControl('emissionColor', label='Color')
        self.endLayout() # end Emission

        self.beginLayout('IDs', collapse=True)
        for i in range(1,9):
            self.addControl('id'+str(i))
        self.endLayout() # end IDs

        self.beginLayout('AOVs', collapse=True)
        self.addControl('lightGroupsIndirect', label='Indirect light groups')
        self.addControl('standardCompatibleAOVs', label='Write standard AOVs only')
        self.addControl('transmitAovs', label='Transmit AOVs')
        self.addControl('aov_diffuse_color', label='Diffuse color')
        self.addControl('aov_direct_diffuse', label='Direct diffuse')
        self.addControl('aov_direct_diffuse_raw', label='Direct diffuse (raw)')
        self.addControl('aov_indirect_diffuse', label='Indirect diffuse')
        self.addControl('aov_indirect_diffuse_raw', label='Indirect diffuse (raw)')
        self.addControl('aov_direct_backlight', label='Direct backlight')
        self.addControl('aov_indirect_backlight', label='Indirect backlight')
        self.addControl('aov_direct_specular', label='Direct specular')
        self.addControl('aov_indirect_specular', label='Indirect specular')
        self.addControl('aov_direct_specular_2', label='Direct specular 2')
        self.addControl('aov_indirect_specular_2', label='Indirect specular 2')
        self.addControl('aov_single_scatter', label='Single scatter')
        self.addControl('aov_sss', label='SSS')
        self.addControl('aov_refraction', label='Refraction')
        self.addControl('aov_emission', label='Emission')
        self.addControl('aov_uv', label='UV')
        self.addControl('aov_depth', label='Depth')
        self.addControl('aov_light_group_1', label='Light group [1]')
        self.addControl('aov_light_group_2', label='Light group [2]')
        self.addControl('aov_light_group_3', label='Light group [3]')
        self.addControl('aov_light_group_4', label='Light group [4]')
        self.addControl('aov_light_group_5', label='Light group [5]')
        self.addControl('aov_light_group_6', label='Light group [6]')
        self.addControl('aov_light_group_7', label='Light group [7]')
        self.addControl('aov_light_group_8', label='Light group [8]')
        self.addControl('aov_id_1', label='ID [1]')
        self.addControl('aov_id_2', label='ID [2]')
        self.addControl('aov_id_3', label='ID [3]')
        self.addControl('aov_id_4', label='ID [4]')
        self.addControl('aov_id_5', label='ID [5]')
        self.addControl('aov_id_6', label='ID [6]')
        self.addControl('aov_id_7', label='ID [7]')
        self.addControl('aov_id_8', label='ID [8]')
        self.endLayout() # end AOVs

        self.addBumpLayout()

        pm.mel.AEdependNodeTemplate(self.nodeName)

        self.addExtraControls()
        self.endScrollLayout()