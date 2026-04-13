import os
import glob
import numpy as np
import nibabel as nib
import pandas as pd
from inspectro_gadget import gadget
from nibabel.affines import apply_affine
import traceback

# =============================================================================
# 1. 경로 설정
# =============================================================================
BASE_DIR = r'D:\RUN'
FIXED_DIR = os.path.join(BASE_DIR, 'Fixed_nii')
INSP_OUTPUT_ROOT = os.path.join(BASE_DIR, 'InSp')

# 폴더 생성
os.makedirs(FIXED_DIR, exist_ok=True)
os.makedirs(INSP_OUTPUT_ROOT, exist_ok=True)

# =============================================================================
# 2. 핵심 처리 함수
# =============================================================================
def process_subject(sub_id):
    print(f"\n{'='*70}")
    print(f"Processing Pipeline for Subject: {sub_id}")
    print(f"{'='*70}")

    # -------------------------------------------------------------------------
    # STEP 1: 마스크 파일 탐색 및 보정 (Fixing)
    # -------------------------------------------------------------------------
    search_pattern = os.path.join(BASE_DIR, f"*{sub_id}*_mask.nii")
    original_files = glob.glob(search_pattern)
    original_files = sorted(list(set(original_files)))
    
    if not original_files:
        print(f"[Skip] No mask files found for {sub_id}.")
        return

    fixed_files_list = []
    roi_labels = []

    for fpath in original_files:
        fname = os.path.basename(fpath)
        
        # 라벨(P번호) 추출
        parts = fname.split('_')
        label = "Unknown"
        for p in parts:
            if p.startswith('P') and p[1:].isdigit():
                label = p
                break
        
        save_path = os.path.join(FIXED_DIR, f"Fixed_{fname}")
        
        try:
            img = nib.load(fpath)
            data = img.get_fdata()
            
            # 0.001 이상 값을 1로 변환 (빈 마스크 체크)
            if np.sum(data > 0.001) == 0:
                print(f"  [Warn] Skipping empty mask: {fname}")
                continue

            new_data = np.zeros(data.shape)
            new_data[data > 0.001] = 1
            new_data = new_data.astype(np.int16)
            
            new_img = nib.Nifti1Image(new_data, img.affine, img.header)
            nib.save(new_img, save_path)
            
            fixed_files_list.append(save_path)
            roi_labels.append(label)
            
        except Exception as e:
            print(f"  [Error] Fix failed for {fname}: {e}")

    if not fixed_files_list:
        print(f"  [Skip] No valid files prepared for {sub_id}.")
        return

    # -------------------------------------------------------------------------
    # STEP 2: Inspectro-Gadget 실행
    # -------------------------------------------------------------------------
    sub_out_dir = os.path.join(INSP_OUTPUT_ROOT, sub_id)
    os.makedirs(sub_out_dir, exist_ok=True)
    
    print(f"  [Gadget] Running analysis for {len(fixed_files_list)} ROIs...")
    
    gadget_result = None
    try:
        gadget_result = gadget.gadget(
            mask_fnames=fixed_files_list,
            mask_labels=roi_labels,
            out_root=sub_out_dir
        )
    except Exception as e:
        print(f"  [Fail] Gadget execution error: {e}")
        return

    # -------------------------------------------------------------------------
    # STEP 3: 유전자 데이터 추출 + 좌표 수동 계산 + 병합
    # -------------------------------------------------------------------------
    print(f"  [Process] Extracting Genes and Calculating MNI Coordinates...")
    
    final_merged_data = [] 
    
    # 유전자 데이터(receptor_data) 가져오기
    if hasattr(gadget_result, 'receptor_data') and isinstance(gadget_result.receptor_data, dict):
        
        for label in roi_labels:
            # A. 유전자 데이터
            expr_df = gadget_result.receptor_data.get(label)
            
            if expr_df is None or expr_df.empty:
                print(f"    [Warn] No gene data for ROI: {label}")
                continue
                
            # B. 해당 ROI 마스크 파일 찾기 (좌표 계산용)
            target_nii_path = None
            for fp in fixed_files_list:
                if label in os.path.basename(fp):
                    target_nii_path = fp
                    break
            
            if target_nii_path is None:
                print(f"    [Error] Mask file missing for {label}")
                continue

            try:
                # C. NIfTI 로드 및 좌표(MNI) 계산
                img = nib.load(target_nii_path)
                data = img.get_fdata()
                affine = img.affine
                
                # 데이터 순서대로 인덱스 추출
                indices = np.where(data > 0.001)
                voxel_indices = np.vstack(indices).T 
                
                # 행 개수 검증
                if len(voxel_indices) != len(expr_df):
                    print(f"    [Mismatch] {label}: Mask({len(voxel_indices)}) != Data({len(expr_df)})")
                    coords_df = pd.DataFrame() 
                else:
                    # 좌표 변환 수행
                    mni_coords = apply_affine(affine, voxel_indices)
                    coords_df = pd.DataFrame(mni_coords, columns=['mni_x', 'mni_y', 'mni_z'])

                # D. 병합 (좌표 + 유전자)
                current_roi_df = pd.concat([coords_df.reset_index(drop=True), 
                                            expr_df.reset_index(drop=True)], axis=1)
                
                # E. 식별자 추가
                current_roi_df.insert(0, 'ROI_Label', label)
                current_roi_df.insert(0, 'Subject_ID', sub_id)
                
                final_merged_data.append(current_roi_df)
                
            except Exception as calc_err:
                print(f"    [Error] Calc failed for {label}: {calc_err}")

    # -------------------------------------------------------------------------
    # STEP 4: 최종 파일 저장
    # -------------------------------------------------------------------------
    if final_merged_data:
        # 모든 ROI 데이터 합치기
        final_df = pd.concat(final_merged_data, ignore_index=True)
        
        # 파일명: Sxxx_InSp_Final_Coords.tsv
        save_path = os.path.join(BASE_DIR, f"{sub_id}_InSp_Final_Coords.tsv")
        final_df.to_csv(save_path, sep='\t', index=False)
        
        print(f"  [SUCCESS] Created: {os.path.basename(save_path)}")
        print(f"  -> Rows: {len(final_df)}, Cols: {len(final_df.columns)}")
        
        has_mni = 'mni_x' in final_df.columns
        print(f"  -> MNI Coordinates Included: {has_mni}")
    else:
        print("  [Fail] No valid data extracted for this subject.")

# =============================================================================
# 3. 전체 실행 루프 (S001 ~ S023)
# =============================================================================
if __name__ == "__main__":
    print("Starting Batch Analysis for S001 - S023...")
    
    for i in range(1, 24):
        try:
            sub_id = f"S{i:03d}"
            # [수정됨] 함수 이름을 정의된 이름(process_subject)과 일치시켰습니다.
            process_subject(sub_id)
        except Exception as e:
            print(f"\n[CRITICAL ERROR] Subject {sub_id} crashed completely.")
            print(traceback.format_exc())
            continue 

    print(f"\n{'='*70}\nAll Tasks Completed.\n{'='*70}")
