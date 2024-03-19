% MSWCMP  1 ���������M���̃E�F�[�u���b�g���g�������k
%
%   MSWCMP �́A�X���b�V���z�[���h�l���v�Z���A�I�������I�v�V�����ɏ]����
%   �E�F�[�u���b�g���g���� 1-D �M���̈��k���s�Ȃ��܂��B
%   OUTPUTS = MSWCMP(OPTION,INPUTS) �́A��ʓI�ȃV���^�b�N�X�ŁAOPTION 
%   �ɑ΂���L���Ȓl�͂��̒ʂ�ł��B
%           'cmp' , 'cmpsig' , 'cmpdec' , 'thr'
%
%   [XC,DECCMP,THRESH] = MSWCMP('cmp',DEC,METH) �܂���  
%   [XC,DECCMP,THRESH] = MSWCMP('cmp',DEC,METH,PARAM) �́A�I���W�i����
%   �����M���s�� X �̈��k�o�[�W���� XC ��Ԃ��܂��B�����ŁA�E�F�[�u���b�g
%   �����̍\���̂� DEC �ł��BXC �́A�E�F�[�u���b�g�W�����X���b�V���z�[���h��
%   ���邱�Ƃœ����܂��BDECCMP �́AXC �Ɋւ���E�F�[�u���b�g������
%   (MDWTDEC ���Q��)�ATHRESH �̓X���b�V���z�[���h�l�̃x�N�g���ł��B
%
%   METH �͈��k��@�̖��O�ŁA�K�v�ȏꍇ�APARAM �͊֘A����p�����[�^�ł�
%   (�ȉ����Q��)�B�o�͈�����I�����邽�߂� 'cmpsig', 'cmpdec', 'thr' ��
%   �g�p���邱�Ƃ��\�ł��B
%       [XC,THRESH] = MSWCMP('cmpsig', ...) �܂���
%       [DECCMP,THRESH] = MSWCMP('cmpdec',...)
%       THRESH = MSWCMP('thr',...) �́A�v�Z�����X���b�V���z�[���h��
%       �Ԃ��܂����A���k�͍s�Ȃ��܂���B
%   
%   �����̍\���̂̓��͈��� DEC �́A4 �̈����Œu���������܂��B
%   DIRDEC, X, WNAME, LEV.
%       [...] = MSWCMP(OPTION,DIRDEC,X,WNAME,LEV,METH) �܂���
%       [...] = MSWCMP(OPTION,DIRDEC,X,WNAME,LEV,METH,PARAM)
%   ���k�����s����A�܂��̓X���b�V���z�[���h���v�Z����O�ɁA�����M���s�� 
%   X �́ADIRDEC �����ɃE�F�[�u���b�g WNAME ���g���ă��x�� LEV �ŕ���
%   ����܂��B
%
%   ����� 3 �̃I�v�V�������͂��g�p�ł��܂��B
%       [...] = MSWCMP(...,S_OR_H) or
%       [...] = MSWCMP(...,S_OR_H,KEEPAPP) or
%       [...] = MSWCMP(...,S_OR_H,KEEPAPP,IDXSIG)
%       - S_OR_H  ('s' �܂��� 'h') �́A�\�t�g�A�܂��̓n�[�h�X���b�V��
%         �z�[���h���Ƃ����Ӗ��ł� (�ڍׂɂ��Ă� MSWTHRESH ���Q��)�B
%       - KEEPAPP (true �܂��� false)�BKEEPAPP �� true �ɓ������ꍇ�A
%         approximation �W���͕ێ�����܂��B
%       - IDXSIG �́A�����M���̃C���f�b�N�X���܂ރx�N�g���A�܂���
%         ������ 'all' �ł��B
%   �f�t�H���g�́A���ꂼ����̂悤�ɂȂ�܂�: 'h', false, 'all'
%
%   �L���Ȉ��k��@ METH �Ɗ֘A����p�����[�^ PARAM �͂��̒ʂ�ł��B
%       'rem_n0'     (0 �ߖT���폜)        
%       'bal_sn'     (�o�����X�X�p�[�X���m����)
%       'sqrtbal_sn' (�o�����X�X�p�[�X���m���� (���))
%   
%       'scarce'     (���k)       
%       'scarcehi'   (�����k) ,  2.5 <= PARAM <= 10
%       'scarceme'   (�����k) ,  1.5 <= PARAM <= 2.5
%       'scarcelo'   (�ሳ�k) ,    1 <= PARAM <= 2
%       PARAM �̓X�p�[�X���p�����[�^�ŁA1 <= PARAM <= 10 �ł���K�v������܂��B
%       ���k�̕��@�ɑ΂��āA����͍s�Ȃ��܂���B
%    
%       'L2_perf'    (�G�l���M�[��)
%       'N0_perf'    (0 �W����)
%       �p�����[�^ PARAM �͕K�v�Ȑ��\��\�������ł� (0 <= PARAM <= 100)�B
%    
%       'glb_thr'    (�O���[�o���ȃX���b�V���z�[���h)
%       �p�����[�^ PARAM �͐��̎����ł��B
%   
%       'man_thr'     (�蓮�̎�@)
%       �p�����[�^ PARAM �́A���̂悤�� NbSIG �s NbLEV ��̍s��A�܂��� 
%       NbSIG �s (NbLEV+1) ��̍s��ł��B
%        - PARAM(i,j) �́Ai �Ԗڂ̐M���ɑ΂��郌�x�� j �� detail �W����
%          �X���b�V���z�[���h (1 <= j <= NbLEV)�B
%        - PARAM(i,NbLEV+1) �́Ai �Ԗڂ̐M���ɑ΂��� approximation �W����
%          �X���b�V���z�[���h (KEEPAPP �� 0 �̏ꍇ)�B
%       �����ŁANbSIG �͐M���̐��ŁANbLEV �͕����̃��x�����ł��B
%
%   �Q�l mdwtdec, mdwtrec, mswthresh, wthresh


%   Copyright 1995-2007 The MathWorks, Inc.
