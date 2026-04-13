%% Gannet Batch Analysis & File Organization
% 1. TSV 결과는 D:\RUN 폴더에 저장
% 2. PDF 결과는 D:\RUN\S0xx 폴더에 정리하여 저장
% 3. Mask 파일은 D:\RUN 폴더에 이름 변경하여 저장

%-----------------------------------------------------------------------
% 초기 설정
%-----------------------------------------------------------------------
Base_Dir = 'D:\Young\MRS_Bialteral_ATL';
Run_Dir = 'D:\RUN';

if ~exist(Run_Dir, 'dir')
    mkdir(Run_Dir);
end

% SPM 및 Gannet 초기화
spm_jobman('initcfg'); 

%-----------------------------------------------------------------------
% Subject Loop (S001 ~ S023)
%-----------------------------------------------------------------------
for sub_idx = 1:23
    % 메모리 및 그래픽 핸들 초기화 (오류 방지 핵심 코드)
    close all; 
    
    % Subject 폴더 경로 설정
    Sub_Name = sprintf('S%03d', sub_idx);
    Data_Dir = fullfile(Base_Dir, Sub_Name);
    
    if ~exist(Data_Dir, 'dir')
        fprintf('Folder not found: %s\n', Data_Dir);
        continue;
    end
    
    fprintf('\nProcessing %s ...\n', Sub_Name);
    cd(Data_Dir);
    
    %-------------------------------------------------------------------
    % 1. 파일 자동 탐색
    %-------------------------------------------------------------------
    P_Files_Struct = dir(fullfile(Data_Dir, 'P*.7'));
    if isempty(P_Files_Struct)
        fprintf('No P-files found in %s\n', Sub_Name);
        continue;
    end
    metabfiles = fullfile(Data_Dir, {P_Files_Struct.name})';
    
    % Anatomy T1 file 찾기
    Anat_Files_Struct = dir(fullfile(Data_Dir, 's*.nii'));
    anat_found = false;
    anat_file = '';
    
    for k = 1:length(Anat_Files_Struct)
        fname = Anat_Files_Struct(k).name;
        if ~contains(fname, 'seg8') && ~startsWith(fname, 'y_') && ~startsWith(fname, 'w')
            anat_file = fullfile(Data_Dir, fname);
            anat_found = true;
            break;
        end
    end
    
    if ~anat_found
        fprintf('No Anatomy file found in %s\n', Sub_Name);
        continue;
    end
    
    anatfiles = repmat({anat_file}, length(metabfiles), 1);
    
    %-------------------------------------------------------------------
    % 2. Gannet Pipeline 구동 (오류 방지를 위해 단계별 close all 추가)
    %-------------------------------------------------------------------
    try
        MRS_struct = GannetLoad(metabfiles);
        close all; % 그림 닫기
        
        MRS_struct = GannetFit(MRS_struct);
        close all;
        
        MRS_struct = GannetCoRegister(MRS_struct, anatfiles);
        close all;
        
        MRS_struct = GannetSegment(MRS_struct);
        close all;
        
    catch ME
        fprintf('Error in Gannet processing for %s: %s\n', Sub_Name, ME.message);
        continue;
    end
    
    %-------------------------------------------------------------------
    % 3. TSV 결과 저장 (D:\RUN\S0xx_Gannet_Results.tsv)
    %-------------------------------------------------------------------
    Filenames = MRS_struct.metabfile';
    
    % [기본] GABA, Glx
    GABA_Cr = MRS_struct.out.vox1.GABA.ConcCr';
    GABA_IU = MRS_struct.out.vox1.GABA.ConcIU';
    Glx_Cr = MRS_struct.out.vox1.Glx.ConcCr';
    Glx_IU = MRS_struct.out.vox1.Glx.ConcIU';
    
    % [추가] NAA, Choline (Cho), Creatine (Cr)
    % 필드가 존재하는지 확인 후 가져오기 (에러 방지)
    if isfield(MRS_struct.out.vox1, 'NAA')
        NAA_Cr = MRS_struct.out.vox1.NAA.ConcCr';
        NAA_IU = MRS_struct.out.vox1.NAA.ConcIU';
    else
        NAA_Cr = nan(length(Filenames), 1);
        NAA_IU = nan(length(Filenames), 1);
    end
    
    if isfield(MRS_struct.out.vox1, 'Cho')
        Cho_Cr = MRS_struct.out.vox1.Cho.ConcCr';
        Cho_IU = MRS_struct.out.vox1.Cho.ConcIU';
    else
        Cho_Cr = nan(length(Filenames), 1);
        Cho_IU = nan(length(Filenames), 1);
    end
    
    % Creatine 자체 농도 (IU)
    if isfield(MRS_struct.out.vox1, 'Cr')
        Cr_IU = MRS_struct.out.vox1.Cr.ConcIU';
    else
        Cr_IU = nan(length(Filenames), 1);
    end

    % 품질 지표 및 조직 분할
    GABA_FitError = MRS_struct.out.vox1.GABA.FitError';
    GABA_FWHM = MRS_struct.out.vox1.GABA.FWHM';
    GM_Fraction = MRS_struct.out.vox1.tissue.fGM';
    WM_Fraction = MRS_struct.out.vox1.tissue.fWM';
    CSF_Fraction = MRS_struct.out.vox1.tissue.fCSF';

    % [통합 테이블 생성]
    ResultsTable = table(Filenames, ...
        GABA_Cr, GABA_IU, ...
        Glx_Cr, Glx_IU, ...
        NAA_Cr, NAA_IU, ...   % 추가됨
        Cho_Cr, Cho_IU, ...   % 추가됨
        Cr_IU, ...            % 추가됨
        GABA_FitError, GABA_FWHM, ...
        GM_Fraction, WM_Fraction, CSF_Fraction);

    % 파일명: S0xx_Gannet_Results_Extended.tsv (기존 파일 덮어쓰지 않게 이름 변경 추천)
    TSV_Name = sprintf('%s_Gannet_Results_Extended.tsv', Sub_Name);
    writetable(ResultsTable, fullfile(Run_Dir, TSV_Name), ...
        'FileType', 'text', 'Delimiter', '\t');
    
    fprintf('  Saved Extended TSV: %s\n', TSV_Name);
    %-------------------------------------------------------------------
    % 4. PDF 결과 이동 (D:\RUN\S0xx 폴더 생성 후 저장) - [요청하신 기능]
    %-------------------------------------------------------------------
    % 저장할 폴더 생성 (예: D:\RUN\S001)
    Sub_PDF_Dir = fullfile(Run_Dir, Sub_Name);
    if ~exist(Sub_PDF_Dir, 'dir')
        mkdir(Sub_PDF_Dir);
    end
    
    % Gannet이 생성하는 4가지 출력 폴더 이름
    Gannet_Folders = {'GannetLoad_output', 'GannetFit_output', 'GannetCoRegister_output', 'GannetSegment_output'};
    
    for gf = 1:length(Gannet_Folders)
        Source_Folder = fullfile(Data_Dir, Gannet_Folders{gf});
        
        if exist(Source_Folder, 'dir')
            % 해당 폴더 내의 모든 PDF 찾기
            Pdf_Files = dir(fullfile(Source_Folder, '*.pdf'));
            
            for p = 1:length(Pdf_Files)
                Src_Pdf = fullfile(Source_Folder, Pdf_Files(p).name);
                Dest_Pdf = fullfile(Sub_PDF_Dir, Pdf_Files(p).name);
                
                % 파일 복사 (원본 유지를 위해 copyfile 사용, 이동하려면 movefile)
                copyfile(Src_Pdf, Dest_Pdf);
            end
        end
    end
    fprintf('Copied PDFs to: %s\n', Sub_PDF_Dir);

    %-------------------------------------------------------------------
    % 5. Inspectro-Gadget용 마스크 생성 및 이동
    %-------------------------------------------------------------------
    def_file = spm_select('FPList', Data_Dir, '^y_.*\.nii$');
    mask_files = spm_select('FPList', Data_Dir, '^P.*_mask\.nii$');
    
    if ~isempty(def_file) && ~isempty(mask_files)
        matlabbatch = {};
        matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(def_file);
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr(mask_files);
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-90 -126 -72; 90 90 108];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 0; 
        
        spm_jobman('run', matlabbatch);
        
        % 생성된 wP*_mask.nii 파일을 RUN 폴더로 이동 및 이름 변경
        Generated_Masks = dir(fullfile(Data_Dir, 'wP*_mask.nii'));
        for m = 1:length(Generated_Masks)
            Src_File = fullfile(Data_Dir, Generated_Masks(m).name);
            Original_Name = Generated_Masks(m).name(2:end); 
            New_Name = sprintf('%s_%s', Sub_Name, Original_Name);
            Dest_File = fullfile(Run_Dir, New_Name);
            
            movefile(Src_File, Dest_File);
        end
    end
end

disp('All processing completed successfully.');
