%GET  WDECTREE �I�u�W�F�N�g�t�B�[���h�̓��e���擾
%
%   [FieldValue1,FieldValue2, ...] = ...
%       GET(T,'FieldName1','FieldName2', ...) 
%   �́AWDECTREE �I�u�W�F�N�g T �Ŏw�肳���t�B�[���h�̓��e��Ԃ��܂��B
%   �I�u�W�F�N�g�A�܂��͍\���̂̃t�B�[���h�ɑ΂��āA�T�u�t�B�[���h�̓��e��
%   �擾���܂� (DTREE/GET ���Q��)�B
%
%   [...] = GET(T) �́AT �̂��ׂẴt�B�[���h�̓��e��Ԃ��܂��B
%
%   'FieldName' �őI���ł���l�͈ȉ��̂Ƃ���ł��B
%   'dtree' - �e�I�u�W�F�N�g
%   'typData' - �f�[�^�̃^�C�v
%   'dimData' - �f�[�^�̎���
%   'WT_Settings' - �E�F�[�u���b�g�ϊ��ݒ�̍\����
%     'typeWT'  - �E�F�[�u���b�g�ϊ��̃^�C�v
%     'wname'   - �E�F�[�u���b�g��
%     'extMode' - DWT �g�����[�h
%     'shift'   - DWT �V�t�g�l
%     'Filters' - �t�B���^�̍\����
%        'Lo_D'    - ���𑤃��[�p�X�t�B���^
%        'Hi_D'    - ���𑤃n�C�p�X�t�B���^
%        'Lo_R'    - �č\�������[�p�X�t�B���^
%        'Hi_R'    - �č\�����n�C�p�X�t�B���^
%
%   �܂��́ADTREE �e�I�u�W�F�N�g�̃t�B�[���h�B
%   'FieldName' �őI���ł���l�ɂ��ẮAhelp dtree/get �Ɠ��͂��Ă��������B
%
%   �Q�l DTREE/READ, DTREE/SET, DTREE/WRITE.


%   Copyright 1995-2008 The MathWorks, Inc.
