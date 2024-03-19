% MSWTHRESH  1 ���������M���̃X���b�V���z�[���h�̎��s
%
%   Y = MSWTHRESH(X,SORH,T) �́A���͍s�� X �́A(SORH = 's' �̏ꍇ) �\�t�g�A
%   �܂��� (SORH = 'h' �̏ꍇ) �n�[�h�� T-�X���b�V���z�[���h��Ԃ��܂��B
%
%   T �́A�P��̒l�AX �Ɠ����T�C�Y�̍s��A�x�N�g���̂����ꂩ�ɂȂ�܂��B
%   �x�N�g���̏ꍇ�A�X���b�V���z�[���h�́A�s�����ɍs�Ȃ��ALT = length(T) 
%   �́ALT => size(X,1) �łȂ���΂Ȃ�܂���B
%
%   Y = MSWTHRESH(X,SORH,T,'c') �́ALT => size(X,2) �ƂȂ�������
%   �X���b�V���z�[���h���s�Ȃ��܂��B
%
%   Y = MSWTHRESH(X,'s',T) �́AY = sign(X).*(abs(X)-T+abs(abs(X)-T))/2 ��
%   �Ԃ��܂��B�\�t�g�X���b�V���z�[���h�ł� 0 �����ɏk������܂��B
%
%   Y = MSWTHRESH(X,'h',T) �́AY = X.*(abs(X)>T) ��Ԃ��܂��B
%   �n�[�h�X���b�V���z�[���h�ł͕s�A���������܂��B
%
%   �Q�l mswden, mswcmp, wthresh, wdencmp, wpdencmp


%   Copyright 1995-2007 The MathWorks, Inc.
