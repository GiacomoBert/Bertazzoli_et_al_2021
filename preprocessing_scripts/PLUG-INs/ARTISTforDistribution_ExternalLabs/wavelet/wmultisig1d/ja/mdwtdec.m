% MDWTDEC  1 ���������M���̃E�F�[�u���b�g����
%
%   DEC = MDWTDEC(DIRDEC,X,LEV,WNAME) �́A�s�� X �ɑ΂��āAWNAME �Ƃ���
%   ���O�̃E�F�[�u���b�g���g���� (DIRDEC = 'r' �̏ꍇ) X �̊e�s�A�܂���
%   (DIRDEC = 'c' �̏ꍇ) �e��̃��x�� LEV �ɂ�����E�F�[�u���b�g������
%   �Ԃ��܂��B
%
%   �o�� DEC �́A���̃t�B�[���h�����\���̂ł��B
%     'dirDec'     : 'r' (�s) �܂��� 'c' (��)
%     'level'      : DWT �����̃��x��
%     'wname'      : �E�F�[�u���b�g��
%     'dwtFilters' : 4 �̃t�B�[���h LoD, HiD, LoR, HiR �����\����
%     'dwtEXTM'    : DWT �g�����[�h (DWTMODE ���Q��)
%     'dwtShift'   : DWT �V�t�g�p�����[�^ (0 �܂��� 1)
%     'dataSize'   : X �̃T�C�Y
%     'ca'         : ���x�� LEV �ɂ����� approximation �W��
%     'cd'         : ���x�� 1 ���烌�x�� LEV �܂ł� detail �W���̃Z���z��
%      �W�� cA �� cD{k} (k = 1 ���� LEV) �͍s��ŁADIRDEC = 'r' (�܂��� 
%      DIRDEC = 'c') �̏ꍇ�A�s���� (�܂��͗����) �Ɋi�[����܂��B
%
%   �E�F�[�u���b�g���̑���ɁA4 �̃t�B���^���g�p���܂��B
%   DEC = MDWTDEC(DIR,X,LEV,LoD,HiD,LoR,HiR)
%
%   MDWTDEC(...,'mode',EXTMODE) �́A�w�肵�� EXTMODE �̊g�����[�h��
%   �E�F�[�u���b�g�������v�Z���܂� (�L���Ȋg�����[�h�ɂ��Ă� DWTMODE 
%   ���Q��)�B
%
%   �Q�l MDWTREC.


%   Copyright 1995-2007 The MathWorks, Inc.
