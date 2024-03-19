% WMULDEN  1 �����E�F�[�u���b�g���ϗʃm�C�Y����
%
%   [X_DEN,NPC,NESTCOV,DEC_DEN,PCA_Params,DEN_Params] = ...
%               WMULDEN(X,LEVEL,WNAME,NPC_APP,NPC_FIN,TPTR,SORH) �܂���
%   [...] = WMULDEN(X,LEVEL,WNAME,'mode',EXTMODE,NPC_APP,...) �́A
%   ���͍s�� X �ɑ΂��ăm�C�Y������̍s�� X_DEN ��Ԃ��܂��B
%   ����́A����m�C�Y�����U�s�񂪁A�Ίp�s��ƂȂ��ϐ��m�C�Y�����@��
%   Approximation �����A�܂��͍č\����̐M���ɑ΂��Ď听�����͂ɂ��
%   �m�C�Y������g�ݍ��킹�����̂ł��B
%
%   ���͍s�� X �́A������Ɋi�[���ꂽ N �̒����̐M�� P ���܂݂܂��B
%   ������ N > P �ł��B
%
%   �E�F�[�u���b�g�����p�����[�^
%   ----------------------------
%   �E�F�[�u���b�g�����́A�������x�� LEVEL �ƃE�F�[�u���b�g WNAME ���g����
%   �s�Ȃ��܂��BEXTMODE �́ADWT �ɑ΂���g�����[�h�ł� (�f�t�H���g�� 
%   DWTMODE �ŕԂ���郂�[�h�ł�)�B
%
%   MDWTDEC ���g���ē����镪�� DEC �����p�\�ȏꍇ�A
%   [...] = WMULDEN(X,LEVEL,WNAME,'mode',EXTMODE,NPC_APP)
%   �̑����  [...] = WMULDEN(DEC,NPC_APP) ���g�p���邱�Ƃ��\�ł��B
%
%   �听���p�����[�^ NPC_APP �� NPC_FIN
%   -----------------------------------
%   ���͂̑I����@ NPC_APP �� NPC_FIN �́A���ꂼ��A�E�F�[�u���b�g�̈��
%   ���x�� LEVEL �ɂ����� Approximation �ƃE�F�[�u���b�g�č\����� PCA ��
%   �΂���听����I��������@���`���܂��B
%
%   NPC_APP (resp. NPC_FIN) �������̏ꍇ�A���x�� LEVEL �ɂ����� approximation 
%   (�܂��̓E�F�[�u���b�g�č\����� PCA) �ɑ΂���c��̎听�������܂݂܂��B
%   NPC_XXX �́A0 <= NPC_XXX <= P �łȂ���΂Ȃ�܂���B
%
%   NPC_APP �܂��� NPC_FIN = 'kais' (�܂��� 'heur') �́AKaiser �̃��[��
%   (�܂��͌o����) ���g���� �c��̎听�����������I�ɑI�����܂��B
%      - Kaiser �̃��[���́A���ׂĂ̌ŗL�l�̕��ς𒴂���ŗL�l�ɑ΂���
%        ������ێ����܂��B
%      - �o�����́A���ׂĂ̌ŗL�l�̘a�� 0.05 �{�𒴂���ŗL�l�ɑ΂���
%        ������ێ����܂��B
%   NPC_APP �܂��� NPC_FIN = 'none' �́ANPC_APP �܂��� NPC_FIN = P �Ɠ����ł��B
%
%   �m�C�Y�����p�����[�^ TPTR, SORH (WDEN �� WBMPEN ���Q��)
%   -------------------------------------------------------
%   �f�t�H���g�l: TPTR = 'sqtwolog' �� SORH = 's'�B
%   TPTR �ɑ΂���L���Ȓl�͂��̒ʂ�ł��B
%       'rigrsure','heursure','sqtwolog','minimaxi'
%       'penalhi','penalme','penallo'
%   SORH �ɑ΂���L���Ȓl�́A's' (�\�t�g) �܂��� 'h' (�n�[�h) �ł��B
%
%   �o��
%   ----
%   X_DEN �́A���͍s�� X �ɑ΂���m�C�Y������̍s��ł��B
%   NPC �́A�I�������c��̎听�����̃x�N�g���ł��B
%   NESTCOV �ŏ��̋����U�̍s�� (MCD) �̐������g���ē����鐄�肳�ꂽ
%   �m�C�Y�̋����U�s��ł��B
%   DEC_DEN �́AX_DEN �̃E�F�[�u���b�g�����ł��B�����̍\���̂̂��ڂ���
%   ���ɂ��Ă� MDWTDEC ���Q�Ƃ��Ă��������B
%   PCA_Params �́A���̂悤�ȍ\���̂ł��B
%       PCA_Params.NEST = {pc_NEST,var_NEST,NESTCOV}
%       PCA_Params.APP  = {pc_APP,var_APP,npc_APP}
%       PCA_Params.FIN  = {pc_FIN,var_FIN,npc_FIN}
%   ������: 
%       - pc_XXX �́A�听���� P�~P �̍s��ł��B
%         ��́A���U���~���ɏ]���Ċi�[���ꂽ���̂ł��B
%       - var_XXX �́A�听���̕��U�x�N�g���ł��B
%       - NESTCOV �́A���x�� 1 �ɂ����� detail �ɑ΂��鋤���U�s��̐���ł��B
%   DEN_Params �́A���̂悤�ȍ\���̂ł��B
%       DEN_Params.thrVAL �́A�e���x���ɑ΂���X���b�V���z�[���h�l��
%       �܂ޒ��� LEVEL �̃x�N�g���ł��B
%       DEN_Params.thrMETH �́A�m�C�Y������@ (TPTR) �̖��O���܂ޕ�����ł��B
%       DEN_Params.thrTYPE �́A�X���b�V���z�[���h (SORH) �̃^�C�v���܂�
%       �����ł��B
%
%   ���ʂȏꍇ
%   ----------
%   [DEC,PCA_Params] = WMULDEN('estimate',DEC,NPC_APP,NPC_FIN) �́A
%   �E�F�[�u���b�g���� DEC �Ǝ听������ PCA_Params ��Ԃ��܂��B
%
%   [X_DEN,NPC,DEC_DEN,PCA_Params] = WMULDEN('execute',DEC,PCA_Params) �܂���
%   [...] = WMULDEN('execute',DEC,PCA_Params,TPTR,SORH) �́A�ȑO�Ɍv�Z���ꂽ
%   �听������ PCA_Params ���g�p���܂��B
%   
%   ���͂̒l DEC �́AX, LEVEL, WNAME �Œu�������邱�Ƃ��\�ł��B


%   Copyright 1995-2007 The MathWorks, Inc.
