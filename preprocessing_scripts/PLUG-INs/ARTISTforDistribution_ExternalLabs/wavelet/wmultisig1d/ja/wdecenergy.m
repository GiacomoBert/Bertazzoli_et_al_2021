% WDECENERGY  1 ���������M���̕����̃G�l���M�[�̕��z
%
%   [E,PEC,PECFS] = WDECENERGY(DEC) �́A�e����M���̃G�l���M�[ (L2-�m����) 
%   ���܂ރx�N�g�� E�A�e�M���̊e�E�F�[�u���b�g���� (approximation �� 
%   detail) �ɑ΂���G�l���M�[�̊������܂ލs�� PEC�A�e�W���ɑ΂���G�l���M�[
%   �̊������܂ލs�� PECFS ���v�Z���܂��B���̂悤�ɂȂ�܂��B
%       - E(i) �́Ai �Ԗڂ̐M���̃G�l���M�[ (L2-�m����) �ł��B
%       - PEC(i,1) �́Ai �Ԗڂ̐M���̃��x�� DEC.level (MAXLEV) ��
%         approximation �ɑ΂���G�l���M�[�̊����ł��B
%       - PEC(i,j) �́Ai �Ԗڂ̐M���̃��x�� (MAXLEV+1-j) �� detail ��
%         �΂���G�l���M�[�̊����ł��B�����ŁAj = 2,...,MAXLEV+1 �ł��B
%       - PECFS(i,j) �́Ai �Ԗڂ̐M���� j �Ԗڂ̌W���ɑ΂���G�l���M�[��
%         �����ł��B
%
%   [E,PEC,PECFS,IDXSORT,LONGS] = WDECENERGY(DEC,'sort') �́A������ 
%   (�s��) �i�[���ꂽ PECFS �ƃC���f�b�N�X�x�N�g�� IDXSORT ��Ԃ��܂��B
%   'sort' �� 'ascend' �Œu��������ƁA�������ʂ������܂��B
%   'sort' �� 'descend' �Œu��������ƁA�~���ɕ��ׂ�ꂽ PECFS �������܂��B
%   LONGS �́A�W���̊e�t�@�~���[�̒������܂ރx�N�g���ł��B
%
%   [...] = WDECENERGY(DEC,OPTSORT,IDXSIG) �́AIDXSIG �x�N�g���ŗ^����ꂽ
%   �C���f�b�N�X�ł���M���ɑ΂���l��Ԃ��܂��BOPTSORT �ɑ΂��ėL����
%   �l�͂��̒ʂ�ł��B
%        'none' , 'sort', 'ascend' , 'descend'.
%
%   �Q�l MDWTDEC, MDWTREC.


%   Copyright 1995-2007 The MathWorks, Inc.
