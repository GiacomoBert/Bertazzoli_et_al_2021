% PLOT   DTREE �I�u�W�F�N�g�̃v���b�g
%
% PLOT(T) �́ADTREE �I�u�W�F�N�g T ���v���b�g���܂��B
% FIG = PLOT(T) �́A�c���[ T ���܂�figure�̃n���h�����o�͂��܂��B
% PLOT(T,FIG)  �́A���Ƀc���[�\�����܂܂��figure FIG �Ƀc���[�\�� T ��
% �v���b�g���܂��B
%
% PLOT �́A�O���t�B�J���ȃc���[�Ǘ����[�e�B���e�B�֐��ł��B
% figure�́A�c���[�\�����܂� GUI �c�[���ł��BNode Level �ł́ADepth_Position �� 
% Index ���ANode Action �ł́ASplit-Merge �� Visualize ��ύX���邱�Ƃ�
% �ł��܂��B�f�t�H���g�l�́ADepth_Position �� Visualize �ł��B
%
% ���݂� Node Action �����s����ɂ́A�m�[�h���N���b�N���܂��B
%
% ���񂩃m�[�h�̕����܂��͑g�ݑւ�������s������ɁA�m�[�h�̊܂܂��
% figure�̃n���h����p���āA�V�����c���[���擾���邱�Ƃ��ł��܂��B
% �ȉ��̓���̃V���^�b�N�X���g��Ȃ���΂Ȃ�܂���:
%       NEWT = PLOT(T,'read',FIG).
%  �����ł́A�ŏ��̈����̓_�~�[�ł��B���̗p�r�Ɋւ���ł���ʓI��
%  �V���^�b�N�X�́A�ȉ��̒ʂ�ł��B:
%       NEWT = PLOT(DUMMY,'READ',FIG);
% �����ŁADUMMY �́ANTREE �I�u�W�F�N�g�ɂ���Đe�ƂȂ�C�ӂ̃I�u�W�F�N�g�ł�
%
% DUMMY �́ANTREE �ɂ���Đe�ƂȂ�I�u�W�F�N�g�ŁA�\�������C�ӂ�
% �I�u�W�F�N�g�����w��ł��܂��B:
%      NEWT = PLOT(ntree,'read',FIG);


%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 12-Mar-96.
%   Last Revision: 13-Feb-2002.
%   Copyright 1995-2004 The MathWorks, Inc.
