rule convert_to_gifti:
    input: join(config['in_freesurfer'],'surf','{hemi}.{surfname}')
    output: 'results/hcp_mmp/sub-{subject}/{hemi}.{surfname}.surf.gii'
    singularity: config['singularity_freesurfer']
    shell: 'mris_convert {input} {output}'

rule convert_to_nifti:
    input: join(config['in_freesurfer'],'mri','{volname}.mgz')
    output: 'results/hcp_mmp/sub-{subject}/{volname}.nii.gz'
    singularity: config['singularity_freesurfer']
    shell: 'mri_convert {input} {output}'



rule get_tkr2scanner:
    input: 
        t1 = 'results/hcp_mmp/sub-{subject}/T1.nii.gz'
    output:
        tkr2scanner = 'results/hcp_mmp/sub-{subject}/tkr2scanner.xfm'
    singularity: config['singularity_freesurfer']
    shell: 'mri_info {input.t1} --tkr2scanner > {output.tkr2scanner}'
     
rule apply_surf_tkr2scanner:
    input: 
        surf = 'results/hcp_mmp/sub-{subject}/{hemi}.{surfname}.surf.gii',
        tkr2scanner = 'results/hcp_mmp/sub-{subject}/tkr2scanner.xfm'
    output: 
        surf = 'results/hcp_mmp/sub-{subject}/{hemi}-scanner.{surfname}.surf.gii'
    threads: 8
    singularity: config['singularity_connectome_workbench']
    shell: 'OMP_NUM_THREADS={threads} wb_command -surface-apply-affine {input.surf} {input.tkr2scanner} {output.surf}'


rule gen_midthickness:
    input:
        white = 'results/hcp_mmp/sub-{subject}/{hemi}.white.surf.gii',
        pial = 'results/hcp_mmp/sub-{subject}/{hemi}.pial.surf.gii'
    output: 
        midthickness = 'results/hcp_mmp/sub-{subject}/{hemi}.midthickness.surf.gii'
    singularity: config['singularity_connectome_workbench']
    threads: 8
    shell: 'OMP_NUM_THREADS={threads} wb_command -surface-average {output.midthickness} -surf {input.white} -surf {input.pial}'
   

rule resample_subj_to_fsaverage_sphere:
    input: 
        surf = 'results/hcp_mmp/sub-{subject}/{hemi}.midthickness.surf.gii',
        current_sphere = 'results/hcp_mmp/sub-{subject}/{hemi}.sphere.reg.surf.gii',
        new_sphere = lambda wildcards: 'resources/standard_mesh_atlases/resample_fsaverage/'
                                        'fs_LR-deformed_to-fsaverage.{H}.sphere.32k_fs_LR.surf.gii'.format(
                                                H = hemi_to_H[wildcards.hemi])
    params:
        method = 'BARYCENTRIC'
    output:
        surf = 'results/hcp_mmp/sub-{subject}/{hemi}.midthickness.32k_fs_LR.surf.gii'
    singularity: config['singularity_connectome_workbench']
    threads: 8
    shell: 'OMP_NUM_THREADS={threads} wb_command -surface-resample {input.surf} {input.current_sphere} {input.new_sphere} {params.method} {output.surf}'


rule resample_labels_to_subj_sphere:
    input:
        label = 'resources/standard_mesh_atlases/{hemi}.hcp-mmp.32k_fs_LR.label.gii',
        current_sphere = lambda wildcards: 'resources/standard_mesh_atlases/resample_fsaverage/'
                                            'fs_LR-deformed_to-fsaverage.{H}.sphere.32k_fs_LR.surf.gii'.format(
                                                H = hemi_to_H[wildcards.hemi]),
        new_sphere = 'results/hcp_mmp/sub-{subject}/{hemi}.sphere.reg.surf.gii',
        current_surf = 'results/hcp_mmp/sub-{subject}/{hemi}.midthickness.32k_fs_LR.surf.gii',
        new_surf = 'results/hcp_mmp/sub-{subject}/{hemi}.midthickness.surf.gii'
    params:
        method = 'ADAP_BARY_AREA'
    output:
        label = 'results/hcp_mmp/sub-{subject}/{hemi}.native.hcp-mmp.label.gii'
    singularity: config['singularity_connectome_workbench']
    threads: 8
    shell: 
        'OMP_NUM_THREADS={threads} wb_command -label-resample {input.label} {input.current_sphere} {input.new_sphere}'
        ' {params.method} {output.label}'
        ' -area-surfs {input.current_surf} {input.new_surf}'



rule map_labels_to_volume_ribbon:
    input: 
        label = 'results/hcp_mmp/sub-{subject}/{hemi}.native.hcp-mmp.label.gii',
        surf = 'results/hcp_mmp/sub-{subject}/{hemi}-scanner.midthickness.surf.gii',
        vol_ref = 'results/hcp_mmp/sub-{subject}/T1.nii.gz',
        white_surf = 'results/hcp_mmp/sub-{subject}/{hemi}-scanner.white.surf.gii',
        pial_surf = 'results/hcp_mmp/sub-{subject}/{hemi}-scanner.pial.surf.gii',
    output:
        label_vol = 'results/hcp_mmp/sub-{subject}/{hemi}.native.hcp-mmp.nii.gz',
    singularity: config['singularity_connectome_workbench']
    threads: 8
    shell:
        'OMP_NUM_THREADS={threads} wb_command -label-to-volume-mapping {input.label} {input.surf} {input.vol_ref} {output.label_vol}'
        ' -ribbon-constrained {input.white_surf} {input.pial_surf}'
        ' -greedy'
     

#currently optional
rule map_labels_to_volume_wmboundary:
    input: 
        label = 'results/hcp_mmp/sub-{subject}/{hemi}.native.hcp-mmp.label.gii',
        surf = 'results/hcp_mmp/sub-{subject}/{hemi}.white.surf.gii',
        vol_ref = 'results/hcp_mmp/sub-{subject}/T1.nii.gz',
        white_surf = 'results/hcp_mmp/sub-{subject}/{hemi}.white.surf.gii',
        pial_surf = 'results/hcp_mmp/sub-{subject}/{hemi}.pial.surf.gii',
    params:
        nearest_vertex = '{wmbdy}' 
    output:
        label_vol = 'results/hcp_mmp/sub-{subject}/{hemi}.native-wmbound{wmbdy}.hcp-mmp.nii.gz',
    singularity: config['singularity_connectome_workbench']
    threads: 8
    shell:
        'OMP_NUM_THREADS={threads} wb_command -label-to-volume-mapping {input.label} {input.surf} {input.vol_ref} {output.label_vol}'
        ' -nearest-vertex {params.nearest_vertex}'
 


