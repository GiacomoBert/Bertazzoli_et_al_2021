% CHGWDECCFS 1 ���������M���̕����W���̕ύX
%
%   DEC = CHGWDECCFS(DEC,'ca',COEFS) �́A���x�� DEC.level �ɂ����� 
%   approximation �W�����s�� COEFS ���Ɋ܂܂��W���ƒu�������܂��B
%   COEFS ���P��̒l V �̏ꍇ�A�W���̂��ׂĂ� V �ɒu���������܂��B
%    
%   DEC = CHGWDECCFS(DEC,'cd',COEFS,LEV) �́A���x�� LEV �ɂ����� 
%   detail �W�����s�� COEFS ���Ɋ܂܂��W���ƒu�������܂��B
%   COEFS ���P��̒l V �̏ꍇ�ALEV �̓��x���̃x�N�g���ŁA������
%   ���x���ɑ����Ă���W���̂��ׂĂ� V �ɒu���������܂��B
%   LEV �� 1 <= LEV <= DEC.level �łȂ���΂Ȃ�܂���B
%
%   DEC = CHGWDECCFS(DEC,'all',CA,CD) �́Aapproximation �W���̂��ׂĂ�
%   detaiil �W���̂��ׂĂ�u�������܂��BCA �͍s��ŁACD �� DEC.level ��
%   �����̃Z���z��łȂ���΂Ȃ�܂���B
%
%   DEC = CHGWDECCFS(...,IDXSIG) �́A�x�N�g�� IDXSIG �ŗ^����ꂽ
%   �C���f�b�N�X�̐M���ɑ΂��ČW����u�������܂��B�����f�[�^���s�� 
%   X �̍s���� (�܂��͗����) �Ɋi�[���ꂽ�ꍇ�AIDXSIG �͊֘A�f�[�^��
%   �s (�܂��͗�) �C���f�b�N�X���܂݂܂��B
%
%   COEFS (or CA, or CD) �� 1 �̐��l (�X�J��) �̏ꍇ�A�֘A�W���̂��ׂĂ�
%   �u�������܂��B�����łȂ���΁ACOEFS (or CA, or CD) �͓K�؂ȃT�C�Y��
%   �s��łȂ���΂Ȃ�܂���B
%
%   ���ۂ̒l V �ɑ΂��āADEC = CHGWDECCFS(DEC,'all',V) �́A�W���̂��ׂĂ�
%   V �ɒu�������܂��B
%
%   �Q�l MDWTDEC, MDWTREC.


%   Copyright 1995-2007 The MathWorks, Inc.
