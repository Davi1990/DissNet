"""
miscelallaneous functions and classes to preprocess raw dwi files
Author: Davide Momi, PhD [momi.davide89@gmail.com], https://twitter.com/davemomi
"""

import glob, os
import nipype.interfaces.mrtrix3 as mrt
import nipype.interfaces.fsl as fsl


class DWI_preproc(object):

    def __init__(self, path, files_name):

        self.path = path
        self.files_name = files_name


    def grab_data(self):
        self.dwi_files=glob.glob(self.path + '/*/' + '*' + self.files_name + '.nii.gz')
        self.bvecs_files=glob.glob(self.path + '/*/' + '*' + 'bvec*')
        self.bvals_files=glob.glob(self.path + '/*/' + '*' + 'bval*')

        return self.dwi_files, self.bvecs_files, self.bvals_files


    def preprocess_dwi_data(self, data, index, acqp):


        if len(data[0]) != len(data[1]):
            raise ValueError('dwi datas do not have the same shape of bvec files')
        if len(data[0]) != len(data[2]):
            raise ValueError('dwi datas do not have the same shape of bval files')
        if len(data[1]) != len(data[2]):
            raise ValueError('bvec files do not have the same shape of bvec files')

        for subj in range(len(data[0])):
            print('Extracting B0 volume for subject',subj)
            self.roi = fsl.ExtractROI(in_file=data[0][subj],
                                      roi_file= os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_nodiff.nii.gz'),
                                      t_min=0, t_size=1)
            self.roi.run()

            print('Converting into .mif for subject',subj)
            self.mrconvert = mrt.MRConvert()
            self.mrconvert.inputs.in_file = data[0][subj]
            self.mrconvert.inputs.grad_fsl = (data[1][subj], data[2][subj])
            self.mrconvert.inputs.out_file = os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_dwi.mif')
            self.mrconvert.run()

            print('Denoising data for subject',subj)
            self.denoise = mrt.DWIDenoise()
            self.denoise.inputs.in_file = os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_dwi.mif')
            self.denoise.inputs.noise = os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_noise.mif')
            self.denoise.inputs.out_file = os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_dwi_denoised.mif')
            self.denoise.run()

            self.denoise_convert = mrt.MRConvert()
            self.denoise_convert.inputs.in_file = os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_dwi_denoised.mif')
            self.denoise_convert.inputs.out_file = os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_dwi_denoised.nii.gz')
            self.denoise_convert.run()

            print('Skull stripping for subject',subj)
            self.mybet = fsl.BET()
            self.mybet.inputs.in_file = os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_nodiff.nii.gz')
            self.mybet.inputs.out_file = os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_denoised_brain.nii.gz')
            self.mybet.inputs.frac = 0.1
            self.mybet.inputs.robust = True
            self.mybet.inputs.mask = True
            self.mybet.run()

            print('Running Eddy for subject',subj)
            self.eddy = Eddy()
            self.eddy.inputs.in_file = os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_dwi_denoised.nii.gz')
            self.eddy.inputs.in_file = os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_denoised_brain_mask.nii.gz')
            self.eddy.inputs.in_acqp = acqp
            self.eddy.inputs.in_bvec = data[1][subj]
            self.eddy.inputs.in_bval = data[2][subj]
            self.eddy.inputs.out_base = os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_dwi_denoised_eddy.nii.gz')
            self.eddy.run()
