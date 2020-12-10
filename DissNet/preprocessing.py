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


    def preprocess_dwi_data(self, data, index, acqp, atlas2use, ResponseSD_algorithm='tournier',
                            fod_algorithm='csd', tract_algorithm='iFOD2',
                            streamlines_number= '10M'):

        '''
        preprocessing of dwi data and connectome extraction

        Parameters
        ----------

        subjects_dir = path to the subjects' folders
        data: tuple |
            a tuple having the path to dwi, bvecs and bvals files. It is obtained
            using the function grab_data()
        index: str |
            Name of text file specifying the relationship between the images in
            --imain and the information in --acqp and --topup. E.g. index.txt
        acqp: str |
            Name of text file with information about the acquisition of the images
            in --imain
        atlas2use: str |
             The input node parcellation image
        ResponseSD_algorithm (optional): str |
             Select the algorithm to be used to complete the script operation;
             Options are: dhollander, fa, manual, msmt_5tt, tax, tournier
             (Default is 'tournier')
        fod_algorithm (optional): str |
             The algorithm to use for FOD estimation. (options are: csd,msmt_csd)
             (Default is 'csd')
        tract_algorithm (optional): str |
            specify the tractography algorithm to use. Valid choices are: FACT,
            iFOD1, iFOD2, Nulldist1, Nulldist2, SD_Stream, Seedtest, Tensor_Det,
            Tensor_Prob (Default is 'iFOD2')
        streamlines_number (optional): str |
            set the desired number of streamlines (Default is '10M')
    '''

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

            print('Running Bias Correction for subject',subj)
            self.bias_correct = mrt.DWIBiasCorrect()
            self.bias_correct.inputs.use_ants = True
            self.bias_correct.inputs.in_file = os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_dwi_denoised_eddy.nii.gz')
            self.bias_correct.inputs.grad_fsl  = (os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_dwi_denoised_eddy.eddy_rotated_bvecs.bvec'),  data[2][subj])
            self.bias_correct.inputs.bias =  os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_bias.mif')
            self.bias_correct.inputs.out_file = os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_dwi_denoised_eddy_unbiased.mif')
            self.bias_correct.run()

            print('Calculating Response function for subject',subj)
            self.resp = mrt.ResponseSD()
            self.resp.inputs.in_file = os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_dwi_denoised_eddy_unbiased.mif')
            self.resp.inputs.algorithm = ResponseSD_algorithm
            self.resp.inputs.grad_fsl = (os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_dwi_denoised_eddy.eddy_rotated_bvecs.bvec'),  data[2][subj])
            self.resp.inputs.wm_file = os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_response.txt')
            self.resp.run()

            print('Estimating FOD for subject',subj)
            self.fod = mrt.EstimateFOD()
            self.fod.inputs.algorithm = fod_algorithm
            self.fod.inputs.in_file = os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_dwi_denoised_eddy_unbiased.mif')
            self.fod.inputs.wm_txt = os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_response.txt')
            self.fod.inputs.mask_file = os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_denoised_brain_mask.nii.gz')
            self.fod.inputs.grad_fsl = os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_response.txt')
            self.fod.run()

            print('Extracting whole brain tract for subject',subj)
            self.tk = mrt.Tractography()
            self.tk.inputs.in_file = os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + 'fods.mif')
            self.tk.inputs.roi_mask = os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_denoised_brain_mask.nii.gz')
            self.tk.inputs.algorithm = tract_algorithm
            self.tk.inputs.seed_image = os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_denoised_brain_mask.nii.gz')
            self.tk.inputs.select = streamlines_number
            self.tk.inputs.out_file = os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_whole_brain_' + streamlines_number + '.tck')
            self.tk.run()

            print('Extracting connectome for subject',subj)
            self.mat = mrt.BuildConnectome()
            self.mat.inputs.in_file = os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_whole_brain_' + streamlines_number + '.tck')
            self.mat.inputs.in_parc = atlas2use
            self.mat.inputs.out_file = os.path.join(os.path.split(data[0][subj])[0] + '/' +  os.path.split(data[0][0])[1].split(".nii.gz")[0] + '_connectome.csv')
            self.mat.run()
