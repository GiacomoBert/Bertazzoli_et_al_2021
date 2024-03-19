% LAURPOLY   �N���X LAURPOLY �ɑ΂���R�X�g���N�^ (Laurent������)
%
% P = LAURPOLY(C,d) �́ALaurent�������I�u�W�F�N�g���o�͂��܂��B
% C �́A�v�f�������� P �̌W���x�N�g���ł���Ad �� P �̒P�����̍ō������ł��B
%
% m ���x�N�g�� C �̒����̏ꍇ�AP �͂���Laurent��������\���܂�:
%     P(z) = C(1)*z^d + C(2)*z^(d-1) + ... + C(m)*z^(d-m+1)
%
% P = LAURPOLY(C,'dmin',d) �́AP �̒P�����̍ō����ł͂Ȃ��A�Œ᎟���w�肵�܂��B
% �Ή�����o�� P �́A����Laurent��������\���܂�:
%     P(z) = C(1)*z^(d+m-1) + ... + C(m-1)*z^(d+1) + C(m)*z^d
%
% P = LAURPOLY(C,'dmax',d) �́AP = LAURPOLY(C,d) �Ɠ����ł��B


%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 19-Mar-2001.
%   Last Revision: 09-Jul-2003.
%   Copyright 1995-2004 The MathWorks, Inc.
