%-----------------------------------------------------------------------
% 초기 설정 (FUS 데이터셋 - 트리구조 오류 수정 버전)
%-----------------------------------------------------------------------
Base_Dir = 'D:\Young\MRS_TUS';  
Run_Dir = 'D:\RUNt';         

if ~exist(Run_Dir, 'dir')
    mkdir(Run_Dir);
end

spm_jobman('initcfg'); 

%-----------------------------------------------------------------------
% Subject Loop (FUS001 ~ FUS024)
%-----------------------------------------------------------------------
for sub_idx = 1:24
    close all; % 메모리 확보
    
    % 폴더명 규칙: FUS001_Sham, FUS002_Sham ...
    Sub_Name = sprintf('FUS%03d_Sham', sub_idx); 
    Data_Dir = fullfile(Base_Dir, Sub_Name);
    
    % 1. 폴더 존재 여부 확인
    if ~exist(Data_Dir, 'dir')
        fprintf('[Skip] Folder not found: %s\n', Sub_Name);
        continue;
    end
    
    fprintf('\n===================================================\n');
    fprintf('Processing %s ...\n', Sub_Name);
    cd(Data_Dir);
    
    %-------------------------------------------------------------------
    % 2. MRS P-files (.7) 찾기
    %-------------------------------------------------------------------
    P_Files_Struct = dir(fullfile(Data_Dir, 'P*.7'));
    if isempty(P_Files_Struct)
        fprintf('[Skip] No P-files found in %s (Folder might be empty)\n', Sub_Name);
        continue; % FUS014_Sham 등 빈 폴더는 여기서 걸러집니다.
    end
    metabfiles = fullfile(Data_Dir, {P_Files_Struct.name})';
    
    %-------------------------------------------------------------------
    % 3. Anatomy T1 file 찾기 (정교한 필터링 + .gz 압축 해제)
    %-------------------------------------------------------------------
    % .nii 파일과 .nii.gz 파일을 모두 찾습니다.
    nii_files = dir(fullfile(Data_Dir, '*.nii'));
    gz_files = dir(fullfile(Data_Dir, '*.nii.gz'));
    All_Anat_Candidates = [nii_files; gz_files];
    
    anat_found = false;
    anat_file = '';
    
    for k = 1:length(All_Anat_Candidates)
        fname = All_Anat_Candidates(k).name;
        
        % [제외 조건] - 결과 파일이나 마스크 파일 등을 제외
        if contains(fname, 'mask') || ...            % 마스크 파일
           startsWith(fname, 'c') || ...             % Segmentation 결과 (c1, c2...)
           startsWith(fname, 'y_') || ...            % Normalization 변형장
           startsWith(fname, 'w') || ...             % Normalized 영상
           startsWith(fname, 'P') || ...             % MRS 마스크 (P로 시작)
           contains(fname, 'seg8') || ...            % Segment info
           All_Anat_Candidates(k).isdir              % 폴더인 경우 제외
       
            continue; 
        end
        
        % [선택 조건] 위 제외 조건을 통과한 파일 중 T1으로 추정되는 것 선택
        % (보통 용량이 크거나 이름에 T1이 들어감)
        
        % 만약 .nii.gz 파일이라면 압축을 해제해야 함
        if endsWith(fname, '.gz')
            fprintf('  Unzipping anatomy file: %s ...\n', fname);
            gunzip(fname);
            fname = fname(1:end-3); % .gz 제거한 이름으로 업데이트
        end
        
        anat_file = fullfile(Data_Dir, fname);
        anat_found = true;
        break; % 첫 번째로 발견된 적절한 파일을 사용
    end
    
    if ~anat_found
        fprintf('[Skip] No valid Anatomy (.nii) file found in %s\n', Sub_Name);
        continue;
    end
    
    fprintf('  Anatomy file used: %s\n', anat_file);
    anatfiles = repmat({anat_file}, length(metabfiles), 1);
    
    %-------------------------------------------------------------------
    % 4. Gannet Pipeline 구동
    %-------------------------------------------------------------------
    try
        % Load
        MRS_struct = GannetLoad(metabfiles);
        close all;
        
        % Fit
        MRS_struct = GannetFit(MRS_struct);
        close all;
        
        % CoRegister (여기서 Mask 생성됨)
        MRS_struct = GannetCoRegister(MRS_struct, anatfiles);
        close all;
        
        % Segment
        MRS_struct = GannetSegment(MRS_struct);
        close all;
        
    catch ME
        fprintf('[Error] Gannet processing failed for %s: %s\n', Sub_Name, ME.message);
        continue;
    end
    
    %-------------------------------------------------------------------
    % 5. TSV 결과 저장
    %-------------------------------------------------------------------
    Filenames = MRS_struct.metabfile';
    
    GABA_Cr = MRS_struct.out.vox1.GABA.ConcCr';
    GABA_IU = MRS_struct.out.vox1.GABA.ConcIU';
    Glx_Cr = MRS_struct.out.vox1.Glx.ConcCr';
    Glx_IU = MRS_struct.out.vox1.Glx.ConcIU';
    
    % [추가] NAA
    if isfield(MRS_struct.out.vox1, 'NAA')
        NAA_Cr = MRS_struct.out.vox1.NAA.ConcCr';
        NAA_IU = MRS_struct.out.vox1.NAA.ConcIU';
    else
        NAA_Cr = nan(length(Filenames), 1);
        NAA_IU = nan(length(Filenames), 1);
    end
    
    % [추가] Choline
    if isfield(MRS_struct.out.vox1, 'Cho')
        Cho_Cr = MRS_struct.out.vox1.Cho.ConcCr';
        Cho_IU = MRS_struct.out.vox1.Cho.ConcIU';
    else
        Cho_Cr = nan(length(Filenames), 1);
        Cho_IU = nan(length(Filenames), 1);
    end
    
    % [추가] Creatine
    if isfield(MRS_struct.out.vox1, 'Cr')
        Cr_IU = MRS_struct.out.vox1.Cr.ConcIU';
    else
        Cr_IU = nan(length(Filenames), 1);
    end

    GABA_FitError = MRS_struct.out.vox1.GABA.FitError';
    GABA_FWHM = MRS_struct.out.vox1.GABA.FWHM';
    GM_Fraction = MRS_struct.out.vox1.tissue.fGM';
    WM_Fraction = MRS_struct.out.vox1.tissue.fWM';
    CSF_Fraction = MRS_struct.out.vox1.tissue.fCSF';

    ResultsTable = table(Filenames, ...
        GABA_Cr, GABA_IU, ...
        Glx_Cr, Glx_IU, ...
        NAA_Cr, NAA_IU, ...
        Cho_Cr, Cho_IU, ...
        Cr_IU, ...
        GABA_FitError, GABA_FWHM, ...
        GM_Fraction, WM_Fraction, CSF_Fraction);

    TSV_Name = sprintf('%s_Gannet_Results_Extended.tsv', Sub_Name);
    writetable(ResultsTable, fullfile(Run_Dir, TSV_Name), ...
        'FileType', 'text', 'Delimiter', '\t');
    
    fprintf('  Saved Extended TSV: %s\n', TSV_Name);
    
    %-------------------------------------------------------------------
    % 6. PDF 결과 이동
    %-------------------------------------------------------------------
    Sub_PDF_Dir = fullfile(Run_Dir, Sub_Name);
    if ~exist(Sub_PDF_Dir, 'dir')
        mkdir(Sub_PDF_Dir);
    end
    
    Gannet_Folders = {'GannetLoad_output', 'GannetFit_output', 'GannetCoRegister_output', 'GannetSegment_output'};
    for gf = 1:length(Gannet_Folders)
        Source_Folder = fullfile(Data_Dir, Gannet_Folders{gf});
        if exist(Source_Folder, 'dir')
            Pdf_Files = dir(fullfile(Source_Folder, '*.pdf'));
            for p = 1:length(Pdf_Files)
                Src_Pdf = fullfile(Source_Folder, Pdf_Files(p).name);
                Dest_Pdf = fullfile(Sub_PDF_Dir, Pdf_Files(p).name);
                copyfile(Src_Pdf, Dest_Pdf);
            end
        end
    end

    %-------------------------------------------------------------------
    % 7. Inspectro-Gadget용 마스크 생성 및 이동
    %-------------------------------------------------------------------
    % Normalization 정보(y_*.nii) 찾기
    def_file_struct = dir(fullfile(Data_Dir, 'y_*.nii'));
    
    if ~isempty(def_file_struct)
        def_file = fullfile(Data_Dir, def_file_struct(1).name);
        mask_files = spm_select('FPList', Data_Dir, '^P.*_mask\.nii$');
        
        if ~isempty(mask_files)
            matlabbatch = {};
            matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(def_file);
            matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr(mask_files);
            matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-90 -126 -72; 90 90 108];
            matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
            matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 0; 
            
            spm_jobman('run', matlabbatch);
            
            % 생성된 wP*_mask.nii 파일을 RUN 폴더로 이동
            Generated_Masks = dir(fullfile(Data_Dir, 'wP*_mask.nii'));
            for m = 1:length(Generated_Masks)
                Src_File = fullfile(Data_Dir, Generated_Masks(m).name);
                Original_Name = Generated_Masks(m).name(2:end); 
                New_Name = sprintf('%s_%s', Sub_Name, Original_Name);
                Dest_File = fullfile(Run_Dir, New_Name);
                
                movefile(Src_File, Dest_File);
                fprintf('  Saved Mask: %s\n', New_Name);
            end
        end
    else
        fprintf('  [Warning] Normalization file (y_*.nii) not found. Skipping mask generation.\n');
    end
end

fprintf('\nAll processing completed.\n');
