% MDWTREC  1 ���������M���̃E�F�[�u���b�g�č\��
%
%   X = MDWTREC(DEC) �́A�E�F�[�u���b�g�����̍\���� DEC ����M����
%   �I���W�i���s���Ԃ��܂� (MDWTDEC ���Q��)�B
%
%   X = MDWTREC(DEC,IDXSIG) �́A�x�N�g�� IDXSIG �ŗ^����ꂽ�C���f�b�N�X��
%   �Ή�����M�����č\�����܂��B
%
%   Y = MDWTREC(DEC,TYPE,LEV) �́A���x�� LEV �ɂ����� detail �܂��� 
%   approximation �W���𒊏o�A�܂��͍č\�����܂� (TYPE �̒l�Ɉˑ�)�B
%       - TYPE �� 'cd' �܂��� 'ca' �̏ꍇ�A���x�� LEV �̌W�������o����܂��B
%       - TYPE �� 'd' �܂��� 'a' �̏ꍇ�A���x�� LEV �̌W�����č\�����ꂽ
%         �M�����o�͂���܂��B
%   LEV �ɑ΂���ő�l�� LEVDEC = DEC.level �ł��B
%
%   TYPE �� 'a' �܂��� 'ca' �̏ꍇ�ALEV �� 0 <= LEV <= LEVDEC �ƂȂ�
%   �悤�Ȑ����łȂ���΂Ȃ�܂���B
%
%   A  = MDWTREC(DEC,'a') �́AA  = MDWTREC(DEC,'a',LEVDEC) �Ɠ����ł��B
%
%   D  = MDWTREC(DEC,'d') �́AX = A + D �ƂȂ�悤�ɁA���ׂĂ� detail 
%   �̘a���܂ލs���Ԃ��܂��B
%
%   CA = MDWTREC(DEC,'ca') �́ACA = MDWTREC(DEC,'ca',LEVDEC) �Ɠ����ł��B
%
%   CD = MDWTREC(DEC,'cd',MODE) �́Adetail �W���̂��ׂĂ��܂ލs���
%   �Ԃ��܂��B
%   CFS = MDWTREC(DEC,'cfs',MODE) �́A�W���̂��ׂĂ��܂ލs���Ԃ��܂��B
%
%   DEC.dirDec �� 'r' (�܂��� 'c') �̏ꍇ�A�A���͍s���� (�܂��͗����) 
%   �ɍs�Ȃ��܂��B
%   MODE = 'descend' (�܂��� 'ascend') �̏ꍇ�A�W���̓��x�� LEVDEC ����
%   1 (�܂��̓��x�� 1 ���烌�x�� LEVDEC) �܂� "�A��" ����܂��B
%   ���͂� MODE ���ȗ����ꂽ�ꍇ�A�f�t�H���g�� MODE = 'descend' �ł��B
%
%   Y = MDWTREC(...,IDXSIG) �́A�x�N�g�� IDXSIG �ŗ^����ꂽ�C���f�b�N�X��
%   �M���ɑ΂��� detail �܂��� approximation �W���𒊏o�A�܂��͍č\��
%   ���܂��B
%
%   �Q�l MDWTDEC.


%   Copyright 1995-2007 The MathWorks, Inc.
