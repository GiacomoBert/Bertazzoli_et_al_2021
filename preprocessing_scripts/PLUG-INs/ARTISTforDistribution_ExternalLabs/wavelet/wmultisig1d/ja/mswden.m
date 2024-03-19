% MSWDEN  1 ���������M���̃E�F�[�u���b�g���g�����m�C�Y����
%
%   MSWDEN �́A�X���b�V���z�[���h���v�Z���A�I�������I�v�V�����ɏ]����
%   �E�F�[�u���b�g���g���� 1-D �M���̃m�C�Y�������s�Ȃ��܂��B
%   OUTPUTS = MSWDEN(OPTION,INPUTS) �́A��ʓI�ȃV���^�b�N�X�ŁAOPTION ��
%   �΂���L���Ȓl�͂��̒ʂ�ł��B
%           'den' , 'densig' , 'dendec' , 'thr'
%
%   [XD,DECDEN,THRESH] = MSWDEN('den',DEC,METH) �܂��� 
%   [XD,DECDEN,THRESH] = MSWDEN('den',DEC,METH,PARAM) �́A
%   �I���W�i���̕����M���s�� X �̃m�C�Y������̍s�� XD ��Ԃ��܂��B
%   �����ŁA�E�F�[�u���b�g�����̍\���̂� DEC �ł��B
%   XD �́A�E�F�[�u���b�g�W�����X���b�V���z�[���h�����邱�Ƃœ����܂��B
%   DECDEN �́AXD �Ɋ֘A����E�F�[�u���b�g������ (MDWTDEC ���Q��)�A
%   THRESH �́A�X���b�V���z�[���h�l�̃x�N�g���ł��B
%
%   METH �́A�m�C�Y������@�̖��O�ŁAPARAM �́A�K�v�ȏꍇ�A�֘A����
%   �p�����[�^�ł� (�ȉ����Q��)�B
%   �o�͈�����I�����邽�߂ɁA'densig' �܂��� 'dendec' �� OPTION ���g�p
%   ���邱�Ƃ��ł��܂��B
%       [XD,THRESH] = MSWDEN('densig', ...) �܂���
%       [DECDEN,THRESH] = MSWDEN('dendec',...)
%
%   �������@�ŁA�X���b�V���z�[���h�l�݂̂����o�����߂ɁA'thr' �� 
%   OPTION ���g�p���邱�Ƃ��ł��܂��B
%   THRESH = MSWDEN('thr',DEC,METH) �܂���
%   THRESH = MSWDEN('thr',DEC,METH,PARAM) �́A�v�Z�����X���b�V���z�[���h��
%   �Ԃ��܂����A�m�C�Y�����͍s�Ȃ��܂���B
%   
%   �����\���̂̓��͈��� DEC �́A4 �̈����Œu���������܂��B
%   DIRDEC, X, WNAME, LEV�B
%       [...] = MSWDEN(OPTION,DIRDEC,X,WNAME,LEV,METH,PARAM).
%   �m�C�Y���������s����A�܂��̓X���b�V���z�[���h���v�Z����O�ɁA
%   �����M�� X �́ADIRDEC �����ɃE�F�[�u���b�g WNAME ���g���ă��x�� 
%   LEV �ŕ�������܂��B
%
%   ����� 3 �̃I�v�V�������͂��g�p�ł��܂��B
%       [...] = MSWDEN(...,S_OR_H) �܂���
%       [...] = MSWDEN(...,S_OR_H,KEEPAPP) �܂���
%       [...] = MSWDEN(...,S_OR_H,KEEPAPP,IDXSIG)
%       - S_or_H  ('s' �܂��� 'h') �́A�\�t�g�A�܂��̓n�[�h�X���b�V��
%         �z�[���h�Ƃ����Ӗ��ł� (�ڍׂɂ��Ă� MSWTHRESH ���Q��)�B
%       - KEEPAPP (true �܂��� false)�BKEEPAPP �� true �ɓ������ꍇ�A
%         approximation �W���͕ێ�����܂��B
%       - IDXSIG �́A�����M���̃C���f�b�N�X���܂ރx�N�g���A�܂���
%         ������ 'all' �ł��B
%   �f�t�H���g�́A���ꂼ����̂悤�ɂȂ�܂�: 'h', false, 'all'.
%
%   �L���ȃm�C�Y������@ METH �Ɗ֘A����p�����[�^ PARAM �͂��̒ʂ�ł��B
%       'rigrsure' Stein Unbiased Risk �̌���
%       'heursure' 'rigrsure' �I�v�V�����̕ό`
%       'sqtwolog' ��ʓI�ȃX���b�V���z�[���h sqrt(2*log(.))
%       'minimaxi' �~�j�}�b�N�X�X���b�V���z�[���h (THSELECT ���Q��)
%       PARAM �́A�ăX�P�[�����O���s���X���b�V���z�[���h�̔{�����`���܂��B
%          'one' �ăX�P�[�����O���s�Ȃ�Ȃ�
%          'sln' ���x�� 1 �̌W���Ɋ�Â����x���m�C�Y�̒P��̐�����g�p����
%                �ăX�P�[�����O���s��
%          'mln' ���x���ˑ��̃��x���m�C�Y�̐�����g�p���čăX�P�[�����O���s��
%
%       'penal'     (�y�i���e�B)       
%       'penalhi'   (�����y�i���e�B)     ,  2.5 <= PARAM <= 10
%       'penalme'   (�����x�̃y�i���e�B) ,  1.5 <= PARAM <= 2.5
%       'penallo'   (�Ⴂ�y�i���e�B)     ,    1 <= PARAM <= 2
%       PARAM �̓X�p�[�X���̃p�����[�^�ŁA1 <= PARAM <= 10 �ł���K�v��
%       ����܂��B�y�i���e�B�̎�@�ɑ΂��āA����͍s�Ȃ��܂���B
%
%       'man_thr'   (�蓮�̎�@)
%       �p�����[�^ PARAM �́A���̂悤�� NbSIG �s NbLEV ��̍s��A
%       �܂��� NbSIG �s (NbLEV+1) ��̍s��ł��B
%         - PARAM(i,j) �́Ai �Ԗڂ̐M���ɑ΂��郌�x�� j �� detail �W����
%          �X���b�V���z�[���h (1 <= j <= NbLEV)�B
%         - PARAM(i,NbLEV+1) �́Ai �Ԗڂ̐M���ɑ΂��� approximation �W����
%          �X���b�V���z�[���h (KEEPAPP �� 0 �̏ꍇ)�B
%
%   �Q�l mdwtdec, mdwtrec, mswthresh, wthresh


%   Copyright 1995-2007 The MathWorks, Inc.
