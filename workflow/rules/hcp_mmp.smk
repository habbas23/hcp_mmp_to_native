rule extract_from_tar:
    input: 
        tar = config['in_freesurfer_tar']
    params:
        out_folder = config['in_freesurfer_root'],
        file_in_tar = 'sub-{subject}/{modality}/{filename}'
    output: 
        filename = join(config['in_freesurfer'],'{modality,surf|mri}','{filename}')
    shell: 'mkdir -p {params.out_folder} && tar -C {params.out_folder} --extract --file={input.tar} {params.file_in_tar}'


def get_gifti_input (wildcards):
    if wildcards.surfname == 'pial': #add .T1 to the name (since pial is a symlink to pial.T1) so can use if extracting from tar
        return join(config['in_freesurfer'],'surf','{hemi}.{surfname}.T1'.format(hemi=H_to_hemi[wildcards.hemi],surfname=wildcards.surfname))
    else:
        return join(config['in_freesurfer'],'surf','{hemi}.{surfname}'.format(hemi=H_to_hemi[wildcards.hemi],surfname=wildcards.surfname))
        
    
rule convert_to_gifti:
    input: get_gifti_input
    output: bids(root='results/hcp_mmp',subject='{subject}',hemi='{hemi}',suffix='{surfname}.surf.gii',space='fsaverage')
    params: 
        license = config['fs_license']
    container: config['singularity_freesurfer']
    log: 'logs/convert_to_gifti/sub-{subject}_{hemi}_{surfname}.log'
    group: 'participant1'
    shell: 'FS_LICENSE={params.license} mris_convert {input} {output} &> {log}'

rule convert_to_nifti:
    input: join(config['in_freesurfer'],'mri','{volname}.mgz')
    output: bids(root='results/hcp_mmp',subject='{subject}',suffix='{volname}.nii.gz')
    params: 
        license = config['fs_license']
    container: config['singularity_freesurfer']
    log: 'logs/convert_to_nifti/sub-{subject}_{volname}.log'
    group: 'participant1'
    shell: 'FS_LICENSE={params.license} mri_convert {input} {output} &> {log}'

rule get_tkr2scanner:
    input: 
        t1 = bids(root='results/hcp_mmp',subject='{subject}',suffix='T1.nii.gz')
    output:
        tkr2scanner = bids(root='results/hcp_mmp',subject='{subject}',suffix='tkr2scanner.xfm')
    params: 
        license = config['fs_license']
    container: config['singularity_freesurfer']
    log: 'logs/get_tkr2scanner/sub-{subject}.log'
    group: 'participant1'
    shell: 'FS_LICENSE={params.license} mri_info {input.t1} --tkr2scanner > {output.tkr2scanner} 2> {log}'
     
rule apply_surf_tkr2scanner:
    input: 
        surf = bids(root='results/hcp_mmp',subject='{subject}',hemi='{hemi}',suffix='{surfname}.surf.gii',space='fsaverage'),
        tkr2scanner = bids(root='results/hcp_mmp',subject='{subject}',suffix='tkr2scanner.xfm')
    output: 
        surf = bids(root='results/hcp_mmp',subject='{subject}',hemi='{hemi}',space='native',suffix='{surfname}.surf.gii')
    threads: 8
    container: config['singularity_connectome_workbench']
    log: 'logs/apply_surf_tkr2scanner/sub-{subject}_{hemi}_{surfname}.log'
    group: 'participant1'
    shell: 'wb_command -surface-apply-affine {input.surf} {input.tkr2scanner} {output.surf} &> {log}'


rule gen_midthickness:
    input:
        white = bids(root='results/hcp_mmp',subject='{subject}',hemi='{hemi}',suffix='white.surf.gii',space='{space}'),
        pial = bids(root='results/hcp_mmp',subject='{subject}',hemi='{hemi}',suffix='pial.surf.gii',space='{space}')
    output: 
        midthickness = bids(root='results/hcp_mmp',subject='{subject}',hemi='{hemi}',suffix='midthickness.surf.gii',space='{space}')
    container: config['singularity_connectome_workbench']
    threads: 8
    log: 'logs/gen_midthickness/sub-{subject}_{hemi}_{space}.log'
    group: 'participant1'
    shell: 'wb_command -surface-average {output.midthickness} -surf {input.white} -surf {input.pial} &> {log}'
   

rule resample_subj_to_fsaverage_sphere:
    input: 
        surf = bids(root='results/hcp_mmp',subject='{subject}',hemi='{hemi}',suffix='midthickness.surf.gii',space='fsaverage'),
        current_sphere = bids(root='results/hcp_mmp',subject='{subject}',hemi='{hemi}',suffix='sphere.reg.surf.gii',space='fsaverage'),
        new_sphere = 'resources/standard_mesh_atlases/resample_fsaverage/'
                                        'fs_LR-deformed_to-fsaverage.{hemi}.sphere.32k_fs_LR.surf.gii'
    params:
        method = 'BARYCENTRIC'
    output:
        surf = bids(root='results/hcp_mmp',subject='{subject}',hemi='{hemi}',space='fsLR',den='32k',suffix='midthickness.surf.gii')
    container: config['singularity_connectome_workbench']
    threads: 8
    log: 'logs/resample_subj_to_fsaverage_sphere/sub-{subject}_{hemi}.log'
    group: 'participant1'
    shell: 'wb_command -surface-resample {input.surf} {input.current_sphere} {input.new_sphere} {params.method} {output.surf} &> {log}'


rule resample_labels_to_subj_sphere:
    input:
        label = lambda wildcards: 'resources/standard_mesh_atlases/{hemi}.hcp-mmp.32k_fs_LR.label.gii'.format(hemi=H_to_hemi[wildcards.hemi]),
        current_sphere = lambda wildcards: 'resources/standard_mesh_atlases/resample_fsaverage/'
                                            'fs_LR-deformed_to-fsaverage.{hemi}.sphere.32k_fs_LR.surf.gii',
        new_sphere = bids(root='results/hcp_mmp',subject='{subject}',hemi='{hemi}',space='fsaverage',suffix='sphere.reg.surf.gii'),
        current_surf = bids(root='results/hcp_mmp',subject='{subject}',hemi='{hemi}',space='fsLR',den='32k',suffix='midthickness.surf.gii'),
        new_surf = bids(root='results/hcp_mmp',subject='{subject}',hemi='{hemi}',space='fsaverage',suffix='midthickness.surf.gii')
    params:
        method = 'ADAP_BARY_AREA'
    output:
        label = bids(root='results/hcp_mmp',subject='{subject}',hemi='{hemi}',label='hcpmmp',space='native',suffix='dseg.label.gii')
    container: config['singularity_connectome_workbench']
    threads: 8
    log: 'logs/resample_labels_to_subj_sphere/sub-{subject}_{hemi}.log'
    group: 'participant1'
    shell: 
        'wb_command -label-resample {input.label} {input.current_sphere} {input.new_sphere}'
        ' {params.method} {output.label}'
        ' -area-surfs {input.current_surf} {input.new_surf} &> {log}'



rule map_labels_to_volume_ribbon:
    input: 
        label = bids(root='results/hcp_mmp',subject='{subject}',hemi='{hemi}',label='hcpmmp',space='native',suffix='dseg.label.gii'),
        surf = bids(root='results/hcp_mmp',subject='{subject}',hemi='{hemi}',space='native',suffix='midthickness.surf.gii'),
        vol_ref = bids(root='results/hcp_mmp',subject='{subject}',suffix='T1.nii.gz'),
        white_surf = bids(root='results/hcp_mmp',subject='{subject}',hemi='{hemi}',space='native',suffix='white.surf.gii'),
        pial_surf = bids(root='results/hcp_mmp',subject='{subject}',hemi='{hemi}',space='native',suffix='pial.surf.gii'),
    output:
        label_vol = bids(root='results/hcp_mmp',subject='{subject}',hemi='{hemi}',label='hcpmmp',space='native',suffix='dseg.nii.gz'),
    container: config['singularity_connectome_workbench']
    threads: 8
    log: 'logs/map_labels_to_volume_ribbon/sub-{subject}_{hemi}.log'
    group: 'participant1'
    shell:
        'wb_command -label-to-volume-mapping {input.label} {input.surf} {input.vol_ref} {output.label_vol}'
        ' -ribbon-constrained {input.white_surf} {input.pial_surf}'
        ' -greedy &> {log}'
     

#currently optional
rule map_labels_to_volume_wmboundary:
    input: 
        label = bids(root='results/hcp_mmp',subject='{subject}',hemi='{hemi}',label='hcpmmp',space='native',suffix='dseg.label.gii'),
        surf = bids(root='results/hcp_mmp',subject='{subject}',hemi='{hemi}',suffix='white.surf.gii',space='native'),
        vol_ref = bids(root='results/hcp_mmp',subject='{subject}',suffix='T1.nii.gz'),
        white_surf = bids(root='results/hcp_mmp',subject='{subject}',hemi='{hemi}',suffix='white.surf.gii',space='native'),
        pial_surf = bids(root='results/hcp_mmp',subject='{subject}',hemi='{hemi}',suffix='pial.surf.gii',space='native'),
    params:
        nearest_vertex = '{wmbdy}' 
    output:
        label_vol = bids(root='results/hcp_mmp',subject='{subject}',hemi='{hemi}',space='native',label='hcpmmp',desc='wmbound{wmbdy}',suffix='dseg.nii.gz'),
    container: config['singularity_connectome_workbench']
    threads: 8
    log: 'logs/map_labels_to_volume_wmboundary/sub-{subject}_{hemi}_wmbound-{wmbdy}.log'
    group: 'participant1'
    shell:
        'wb_command -label-to-volume-mapping {input.label} {input.surf} {input.vol_ref} {output.label_vol}'
        ' -nearest-vertex {params.nearest_vertex} &> {log}'
 


