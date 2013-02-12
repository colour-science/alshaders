import pymel.core as pm
import maya.cmds as cmds
from alShaders import alShadersTemplate

class AEalPhotometricTemplate(alShadersTemplate):
    
    def filenameEdit(self, data):
        attr = self.nodeAttr('filename')
        cmds.setAttr(attr, data)

    def fileDialog(self, *args):
        fn = cmds.fileDialog2(fileFilter='*.ies *.IES', dialogStyle=2, cap='Load IES file', okc='Load', fm=1)
        if fn is not None and len(fn):
            filenameEdit(fn[0])
            cmds.textFieldButtonGrp('filenameGrp', edit=True, text=fn[0])

    def filenameNew(self, nodeName):
        path = cmds.textFieldButtonGrp('filenameGrp', label='IES file', changeCommand=self.filenameEdit, width=300)
        cmds.textFieldButtonGrp(path, edit=True, text=cmds.getAttr(nodeName))
        cmds.textFieldButtonGrp(path, edit=True, buttonLabel='...', buttonCommand=self.fileDialog)

    def filenameReplace(self, nodeName):
        cmds.textFieldButtonGrp('filenameGrp', edit=True, text=cmds.getAttr(nodeName))

    def setup(self):
        self.addSwatch()
        self.beginScrollLayout()

        self.addCustom('filename', self.filenameNew, self.filenameReplace)

        pm.mel.AEdependNodeTemplate(self.nodeName)

        self.addExtraControls()
        self.endScrollLayout()