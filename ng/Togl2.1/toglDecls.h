#ifndef ToglDecls_H
#  define ToglDecls_H

/* 
 * Togl - a Tk OpenGL widget
 *
 * Copyright (C) 1996-2002  Brian Paul and Ben Bederson
 * Copyright (C) 2005-2009  Greg Couch
 * See the LICENSE file for copyright details.
 */

/* !BEGIN!: Do not edit below this line. */

/*
 * Exported function declarations:
 */

#ifndef Togl_Init_TCL_DECLARED
#define Togl_Init_TCL_DECLARED
/* 0 */
EXTERN int		Togl_Init(Tcl_Interp *interp);
#endif
#ifndef Togl_MakeCurrent_TCL_DECLARED
#define Togl_MakeCurrent_TCL_DECLARED
/* 1 */
EXTERN void		Togl_MakeCurrent(const Togl *togl);
#endif
#ifndef Togl_PostRedisplay_TCL_DECLARED
#define Togl_PostRedisplay_TCL_DECLARED
/* 2 */
EXTERN void		Togl_PostRedisplay(Togl *togl);
#endif
#ifndef Togl_SwapBuffers_TCL_DECLARED
#define Togl_SwapBuffers_TCL_DECLARED
/* 3 */
EXTERN void		Togl_SwapBuffers(const Togl *togl);
#endif
#ifndef Togl_Ident_TCL_DECLARED
#define Togl_Ident_TCL_DECLARED
/* 4 */
EXTERN const char *	Togl_Ident(const Togl *togl);
#endif
#ifndef Togl_Width_TCL_DECLARED
#define Togl_Width_TCL_DECLARED
/* 5 */
EXTERN int		Togl_Width(const Togl *togl);
#endif
#ifndef Togl_Height_TCL_DECLARED
#define Togl_Height_TCL_DECLARED
/* 6 */
EXTERN int		Togl_Height(const Togl *togl);
#endif
#ifndef Togl_Interp_TCL_DECLARED
#define Togl_Interp_TCL_DECLARED
/* 7 */
EXTERN Tcl_Interp *	Togl_Interp(const Togl *togl);
#endif
#ifndef Togl_TkWin_TCL_DECLARED
#define Togl_TkWin_TCL_DECLARED
/* 8 */
EXTERN Tk_Window	Togl_TkWin(const Togl *togl);
#endif
#ifndef Togl_CommandName_TCL_DECLARED
#define Togl_CommandName_TCL_DECLARED
/* 9 */
EXTERN const char *	Togl_CommandName(const Togl *togl);
#endif
#ifndef Togl_AllocColor_TCL_DECLARED
#define Togl_AllocColor_TCL_DECLARED
/* 10 */
EXTERN unsigned long	Togl_AllocColor(const Togl *togl, float red,
				float green, float blue);
#endif
#ifndef Togl_FreeColor_TCL_DECLARED
#define Togl_FreeColor_TCL_DECLARED
/* 11 */
EXTERN void		Togl_FreeColor(const Togl *togl, unsigned long index);
#endif
#ifndef Togl_SetColor_TCL_DECLARED
#define Togl_SetColor_TCL_DECLARED
/* 12 */
EXTERN void		Togl_SetColor(const Togl *togl, unsigned long index,
				float red, float green, float blue);
#endif
#ifndef Togl_LoadBitmapFont_TCL_DECLARED
#define Togl_LoadBitmapFont_TCL_DECLARED
/* 13 */
EXTERN Tcl_Obj *	Togl_LoadBitmapFont(const Togl *togl,
				const char *fontname);
#endif
#ifndef Togl_UnloadBitmapFont_TCL_DECLARED
#define Togl_UnloadBitmapFont_TCL_DECLARED
/* 14 */
EXTERN int		Togl_UnloadBitmapFont(const Togl *togl,
				Tcl_Obj *toglfont);
#endif
#ifndef Togl_UseLayer_TCL_DECLARED
#define Togl_UseLayer_TCL_DECLARED
/* 15 */
EXTERN void		Togl_UseLayer(Togl *togl, int layer);
#endif
#ifndef Togl_ShowOverlay_TCL_DECLARED
#define Togl_ShowOverlay_TCL_DECLARED
/* 16 */
EXTERN void		Togl_ShowOverlay(Togl *togl);
#endif
#ifndef Togl_HideOverlay_TCL_DECLARED
#define Togl_HideOverlay_TCL_DECLARED
/* 17 */
EXTERN void		Togl_HideOverlay(Togl *togl);
#endif
#ifndef Togl_PostOverlayRedisplay_TCL_DECLARED
#define Togl_PostOverlayRedisplay_TCL_DECLARED
/* 18 */
EXTERN void		Togl_PostOverlayRedisplay(Togl *togl);
#endif
#ifndef Togl_ExistsOverlay_TCL_DECLARED
#define Togl_ExistsOverlay_TCL_DECLARED
/* 19 */
EXTERN int		Togl_ExistsOverlay(const Togl *togl);
#endif
#ifndef Togl_GetOverlayTransparentValue_TCL_DECLARED
#define Togl_GetOverlayTransparentValue_TCL_DECLARED
/* 20 */
EXTERN int		Togl_GetOverlayTransparentValue(const Togl *togl);
#endif
#ifndef Togl_IsMappedOverlay_TCL_DECLARED
#define Togl_IsMappedOverlay_TCL_DECLARED
/* 21 */
EXTERN int		Togl_IsMappedOverlay(const Togl *togl);
#endif
#ifndef Togl_AllocColorOverlay_TCL_DECLARED
#define Togl_AllocColorOverlay_TCL_DECLARED
/* 22 */
EXTERN unsigned long	Togl_AllocColorOverlay(const Togl *togl, float red,
				float green, float blue);
#endif
#ifndef Togl_FreeColorOverlay_TCL_DECLARED
#define Togl_FreeColorOverlay_TCL_DECLARED
/* 23 */
EXTERN void		Togl_FreeColorOverlay(const Togl *togl,
				unsigned long index);
#endif
#ifndef Togl_GetClientData_TCL_DECLARED
#define Togl_GetClientData_TCL_DECLARED
/* 24 */
EXTERN ClientData	Togl_GetClientData(const Togl *togl);
#endif
#ifndef Togl_SetClientData_TCL_DECLARED
#define Togl_SetClientData_TCL_DECLARED
/* 25 */
EXTERN void		Togl_SetClientData(Togl *togl, ClientData clientData);
#endif
#ifndef Togl_DrawBuffer_TCL_DECLARED
#define Togl_DrawBuffer_TCL_DECLARED
/* 26 */
EXTERN void		Togl_DrawBuffer(Togl *togl, GLenum mode);
#endif
#ifndef Togl_Clear_TCL_DECLARED
#define Togl_Clear_TCL_DECLARED
/* 27 */
EXTERN void		Togl_Clear(const Togl *togl, GLbitfield mask);
#endif
#ifndef Togl_Frustum_TCL_DECLARED
#define Togl_Frustum_TCL_DECLARED
/* 28 */
EXTERN void		Togl_Frustum(const Togl *togl, GLdouble left,
				GLdouble right, GLdouble bottom,
				GLdouble top, GLdouble near, GLdouble far);
#endif
#ifndef Togl_GetToglFromObj_TCL_DECLARED
#define Togl_GetToglFromObj_TCL_DECLARED
/* 29 */
EXTERN int		Togl_GetToglFromObj(Tcl_Interp *interp, Tcl_Obj *obj,
				Togl **toglPtr);
#endif
#ifndef Togl_TakePhoto_TCL_DECLARED
#define Togl_TakePhoto_TCL_DECLARED
/* 30 */
EXTERN int		Togl_TakePhoto(Togl *togl, Tk_PhotoHandle photo);
#endif
#ifndef Togl_GetProcAddr_TCL_DECLARED
#define Togl_GetProcAddr_TCL_DECLARED
/* 31 */
EXTERN Togl_FuncPtr	Togl_GetProcAddr(const char *funcname);
#endif
#ifndef Togl_GetToglFromName_TCL_DECLARED
#define Togl_GetToglFromName_TCL_DECLARED
/* 32 */
EXTERN int		Togl_GetToglFromName(Tcl_Interp *interp,
				const char *cmdName, Togl **toglPtr);
#endif
#ifndef Togl_SwapInterval_TCL_DECLARED
#define Togl_SwapInterval_TCL_DECLARED
/* 33 */
EXTERN Bool		Togl_SwapInterval(const Togl *togl, int interval);
#endif
#ifndef Togl_Ortho_TCL_DECLARED
#define Togl_Ortho_TCL_DECLARED
/* 34 */
EXTERN void		Togl_Ortho(const Togl *togl, GLdouble left,
				GLdouble right, GLdouble bottom,
				GLdouble top, GLdouble near, GLdouble far);
#endif
#ifndef Togl_NumEyes_TCL_DECLARED
#define Togl_NumEyes_TCL_DECLARED
/* 35 */
EXTERN int		Togl_NumEyes(const Togl *togl);
#endif
#ifndef Togl_ContextTag_TCL_DECLARED
#define Togl_ContextTag_TCL_DECLARED
/* 36 */
EXTERN int		Togl_ContextTag(const Togl *togl);
#endif
#ifndef Togl_UpdatePending_TCL_DECLARED
#define Togl_UpdatePending_TCL_DECLARED
/* 37 */
EXTERN Bool		Togl_UpdatePending(const Togl *togl);
#endif
#ifndef Togl_WriteObj_TCL_DECLARED
#define Togl_WriteObj_TCL_DECLARED
/* 38 */
EXTERN int		Togl_WriteObj(const Togl *togl,
				const Tcl_Obj *toglfont, Tcl_Obj *obj);
#endif
#ifndef Togl_WriteChars_TCL_DECLARED
#define Togl_WriteChars_TCL_DECLARED
/* 39 */
EXTERN int		Togl_WriteChars(const Togl *togl,
				const Tcl_Obj *toglfont, const char *str,
				int len);
#endif
#ifndef Togl_HasRGBA_TCL_DECLARED
#define Togl_HasRGBA_TCL_DECLARED
/* 40 */
EXTERN Bool		Togl_HasRGBA(const Togl *togl);
#endif
#ifndef Togl_IsDoubleBuffered_TCL_DECLARED
#define Togl_IsDoubleBuffered_TCL_DECLARED
/* 41 */
EXTERN Bool		Togl_IsDoubleBuffered(const Togl *togl);
#endif
#ifndef Togl_HasDepthBuffer_TCL_DECLARED
#define Togl_HasDepthBuffer_TCL_DECLARED
/* 42 */
EXTERN Bool		Togl_HasDepthBuffer(const Togl *togl);
#endif
#ifndef Togl_HasAccumulationBuffer_TCL_DECLARED
#define Togl_HasAccumulationBuffer_TCL_DECLARED
/* 43 */
EXTERN Bool		Togl_HasAccumulationBuffer(const Togl *togl);
#endif
#ifndef Togl_HasDestinationAlpha_TCL_DECLARED
#define Togl_HasDestinationAlpha_TCL_DECLARED
/* 44 */
EXTERN Bool		Togl_HasDestinationAlpha(const Togl *togl);
#endif
#ifndef Togl_HasStencilBuffer_TCL_DECLARED
#define Togl_HasStencilBuffer_TCL_DECLARED
/* 45 */
EXTERN Bool		Togl_HasStencilBuffer(const Togl *togl);
#endif
#ifndef Togl_StereoMode_TCL_DECLARED
#define Togl_StereoMode_TCL_DECLARED
/* 46 */
EXTERN int		Togl_StereoMode(const Togl *togl);
#endif
#ifndef Togl_HasMultisample_TCL_DECLARED
#define Togl_HasMultisample_TCL_DECLARED
/* 47 */
EXTERN Bool		Togl_HasMultisample(const Togl *togl);
#endif
#ifndef Togl_CopyContext_TCL_DECLARED
#define Togl_CopyContext_TCL_DECLARED
/* 48 */
EXTERN int		Togl_CopyContext(const Togl *from, const Togl *to,
				unsigned int mask);
#endif

typedef struct ToglStubs {
    int magic;
    const struct ToglStubHooks *hooks;

    int (*togl_Init) (Tcl_Interp *interp); /* 0 */
    void (*togl_MakeCurrent) (const Togl *togl); /* 1 */
    void (*togl_PostRedisplay) (Togl *togl); /* 2 */
    void (*togl_SwapBuffers) (const Togl *togl); /* 3 */
    const char * (*togl_Ident) (const Togl *togl); /* 4 */
    int (*togl_Width) (const Togl *togl); /* 5 */
    int (*togl_Height) (const Togl *togl); /* 6 */
    Tcl_Interp * (*togl_Interp) (const Togl *togl); /* 7 */
    Tk_Window (*togl_TkWin) (const Togl *togl); /* 8 */
    const char * (*togl_CommandName) (const Togl *togl); /* 9 */
    unsigned long (*togl_AllocColor) (const Togl *togl, float red, float green, float blue); /* 10 */
    void (*togl_FreeColor) (const Togl *togl, unsigned long index); /* 11 */
    void (*togl_SetColor) (const Togl *togl, unsigned long index, float red, float green, float blue); /* 12 */
    Tcl_Obj * (*togl_LoadBitmapFont) (const Togl *togl, const char *fontname); /* 13 */
    int (*togl_UnloadBitmapFont) (const Togl *togl, Tcl_Obj *toglfont); /* 14 */
    void (*togl_UseLayer) (Togl *togl, int layer); /* 15 */
    void (*togl_ShowOverlay) (Togl *togl); /* 16 */
    void (*togl_HideOverlay) (Togl *togl); /* 17 */
    void (*togl_PostOverlayRedisplay) (Togl *togl); /* 18 */
    int (*togl_ExistsOverlay) (const Togl *togl); /* 19 */
    int (*togl_GetOverlayTransparentValue) (const Togl *togl); /* 20 */
    int (*togl_IsMappedOverlay) (const Togl *togl); /* 21 */
    unsigned long (*togl_AllocColorOverlay) (const Togl *togl, float red, float green, float blue); /* 22 */
    void (*togl_FreeColorOverlay) (const Togl *togl, unsigned long index); /* 23 */
    ClientData (*togl_GetClientData) (const Togl *togl); /* 24 */
    void (*togl_SetClientData) (Togl *togl, ClientData clientData); /* 25 */
    void (*togl_DrawBuffer) (Togl *togl, GLenum mode); /* 26 */
    void (*togl_Clear) (const Togl *togl, GLbitfield mask); /* 27 */
    void (*togl_Frustum) (const Togl *togl, GLdouble left, GLdouble right, GLdouble bottom, GLdouble top, GLdouble near, GLdouble far); /* 28 */
    int (*togl_GetToglFromObj) (Tcl_Interp *interp, Tcl_Obj *obj, Togl **toglPtr); /* 29 */
    int (*togl_TakePhoto) (Togl *togl, Tk_PhotoHandle photo); /* 30 */
    Togl_FuncPtr (*togl_GetProcAddr) (const char *funcname); /* 31 */
    int (*togl_GetToglFromName) (Tcl_Interp *interp, const char *cmdName, Togl **toglPtr); /* 32 */
    Bool (*togl_SwapInterval) (const Togl *togl, int interval); /* 33 */
    void (*togl_Ortho) (const Togl *togl, GLdouble left, GLdouble right, GLdouble bottom, GLdouble top, GLdouble near, GLdouble far); /* 34 */
    int (*togl_NumEyes) (const Togl *togl); /* 35 */
    int (*togl_ContextTag) (const Togl *togl); /* 36 */
    Bool (*togl_UpdatePending) (const Togl *togl); /* 37 */
    int (*togl_WriteObj) (const Togl *togl, const Tcl_Obj *toglfont, Tcl_Obj *obj); /* 38 */
    int (*togl_WriteChars) (const Togl *togl, const Tcl_Obj *toglfont, const char *str, int len); /* 39 */
    Bool (*togl_HasRGBA) (const Togl *togl); /* 40 */
    Bool (*togl_IsDoubleBuffered) (const Togl *togl); /* 41 */
    Bool (*togl_HasDepthBuffer) (const Togl *togl); /* 42 */
    Bool (*togl_HasAccumulationBuffer) (const Togl *togl); /* 43 */
    Bool (*togl_HasDestinationAlpha) (const Togl *togl); /* 44 */
    Bool (*togl_HasStencilBuffer) (const Togl *togl); /* 45 */
    int (*togl_StereoMode) (const Togl *togl); /* 46 */
    Bool (*togl_HasMultisample) (const Togl *togl); /* 47 */
    int (*togl_CopyContext) (const Togl *from, const Togl *to, unsigned int mask); /* 48 */
} ToglStubs;

#if defined(USE_TOGL_STUBS) && !defined(USE_TOGL_STUB_PROCS)
extern const ToglStubs *toglStubsPtr;
#endif /* defined(USE_TOGL_STUBS) && !defined(USE_TOGL_STUB_PROCS) */

#if defined(USE_TOGL_STUBS) && !defined(USE_TOGL_STUB_PROCS)

/*
 * Inline function declarations:
 */

#ifndef Togl_Init
#define Togl_Init \
	(toglStubsPtr->togl_Init) /* 0 */
#endif
#ifndef Togl_MakeCurrent
#define Togl_MakeCurrent \
	(toglStubsPtr->togl_MakeCurrent) /* 1 */
#endif
#ifndef Togl_PostRedisplay
#define Togl_PostRedisplay \
	(toglStubsPtr->togl_PostRedisplay) /* 2 */
#endif
#ifndef Togl_SwapBuffers
#define Togl_SwapBuffers \
	(toglStubsPtr->togl_SwapBuffers) /* 3 */
#endif
#ifndef Togl_Ident
#define Togl_Ident \
	(toglStubsPtr->togl_Ident) /* 4 */
#endif
#ifndef Togl_Width
#define Togl_Width \
	(toglStubsPtr->togl_Width) /* 5 */
#endif
#ifndef Togl_Height
#define Togl_Height \
	(toglStubsPtr->togl_Height) /* 6 */
#endif
#ifndef Togl_Interp
#define Togl_Interp \
	(toglStubsPtr->togl_Interp) /* 7 */
#endif
#ifndef Togl_TkWin
#define Togl_TkWin \
	(toglStubsPtr->togl_TkWin) /* 8 */
#endif
#ifndef Togl_CommandName
#define Togl_CommandName \
	(toglStubsPtr->togl_CommandName) /* 9 */
#endif
#ifndef Togl_AllocColor
#define Togl_AllocColor \
	(toglStubsPtr->togl_AllocColor) /* 10 */
#endif
#ifndef Togl_FreeColor
#define Togl_FreeColor \
	(toglStubsPtr->togl_FreeColor) /* 11 */
#endif
#ifndef Togl_SetColor
#define Togl_SetColor \
	(toglStubsPtr->togl_SetColor) /* 12 */
#endif
#ifndef Togl_LoadBitmapFont
#define Togl_LoadBitmapFont \
	(toglStubsPtr->togl_LoadBitmapFont) /* 13 */
#endif
#ifndef Togl_UnloadBitmapFont
#define Togl_UnloadBitmapFont \
	(toglStubsPtr->togl_UnloadBitmapFont) /* 14 */
#endif
#ifndef Togl_UseLayer
#define Togl_UseLayer \
	(toglStubsPtr->togl_UseLayer) /* 15 */
#endif
#ifndef Togl_ShowOverlay
#define Togl_ShowOverlay \
	(toglStubsPtr->togl_ShowOverlay) /* 16 */
#endif
#ifndef Togl_HideOverlay
#define Togl_HideOverlay \
	(toglStubsPtr->togl_HideOverlay) /* 17 */
#endif
#ifndef Togl_PostOverlayRedisplay
#define Togl_PostOverlayRedisplay \
	(toglStubsPtr->togl_PostOverlayRedisplay) /* 18 */
#endif
#ifndef Togl_ExistsOverlay
#define Togl_ExistsOverlay \
	(toglStubsPtr->togl_ExistsOverlay) /* 19 */
#endif
#ifndef Togl_GetOverlayTransparentValue
#define Togl_GetOverlayTransparentValue \
	(toglStubsPtr->togl_GetOverlayTransparentValue) /* 20 */
#endif
#ifndef Togl_IsMappedOverlay
#define Togl_IsMappedOverlay \
	(toglStubsPtr->togl_IsMappedOverlay) /* 21 */
#endif
#ifndef Togl_AllocColorOverlay
#define Togl_AllocColorOverlay \
	(toglStubsPtr->togl_AllocColorOverlay) /* 22 */
#endif
#ifndef Togl_FreeColorOverlay
#define Togl_FreeColorOverlay \
	(toglStubsPtr->togl_FreeColorOverlay) /* 23 */
#endif
#ifndef Togl_GetClientData
#define Togl_GetClientData \
	(toglStubsPtr->togl_GetClientData) /* 24 */
#endif
#ifndef Togl_SetClientData
#define Togl_SetClientData \
	(toglStubsPtr->togl_SetClientData) /* 25 */
#endif
#ifndef Togl_DrawBuffer
#define Togl_DrawBuffer \
	(toglStubsPtr->togl_DrawBuffer) /* 26 */
#endif
#ifndef Togl_Clear
#define Togl_Clear \
	(toglStubsPtr->togl_Clear) /* 27 */
#endif
#ifndef Togl_Frustum
#define Togl_Frustum \
	(toglStubsPtr->togl_Frustum) /* 28 */
#endif
#ifndef Togl_GetToglFromObj
#define Togl_GetToglFromObj \
	(toglStubsPtr->togl_GetToglFromObj) /* 29 */
#endif
#ifndef Togl_TakePhoto
#define Togl_TakePhoto \
	(toglStubsPtr->togl_TakePhoto) /* 30 */
#endif
#ifndef Togl_GetProcAddr
#define Togl_GetProcAddr \
	(toglStubsPtr->togl_GetProcAddr) /* 31 */
#endif
#ifndef Togl_GetToglFromName
#define Togl_GetToglFromName \
	(toglStubsPtr->togl_GetToglFromName) /* 32 */
#endif
#ifndef Togl_SwapInterval
#define Togl_SwapInterval \
	(toglStubsPtr->togl_SwapInterval) /* 33 */
#endif
#ifndef Togl_Ortho
#define Togl_Ortho \
	(toglStubsPtr->togl_Ortho) /* 34 */
#endif
#ifndef Togl_NumEyes
#define Togl_NumEyes \
	(toglStubsPtr->togl_NumEyes) /* 35 */
#endif
#ifndef Togl_ContextTag
#define Togl_ContextTag \
	(toglStubsPtr->togl_ContextTag) /* 36 */
#endif
#ifndef Togl_UpdatePending
#define Togl_UpdatePending \
	(toglStubsPtr->togl_UpdatePending) /* 37 */
#endif
#ifndef Togl_WriteObj
#define Togl_WriteObj \
	(toglStubsPtr->togl_WriteObj) /* 38 */
#endif
#ifndef Togl_WriteChars
#define Togl_WriteChars \
	(toglStubsPtr->togl_WriteChars) /* 39 */
#endif
#ifndef Togl_HasRGBA
#define Togl_HasRGBA \
	(toglStubsPtr->togl_HasRGBA) /* 40 */
#endif
#ifndef Togl_IsDoubleBuffered
#define Togl_IsDoubleBuffered \
	(toglStubsPtr->togl_IsDoubleBuffered) /* 41 */
#endif
#ifndef Togl_HasDepthBuffer
#define Togl_HasDepthBuffer \
	(toglStubsPtr->togl_HasDepthBuffer) /* 42 */
#endif
#ifndef Togl_HasAccumulationBuffer
#define Togl_HasAccumulationBuffer \
	(toglStubsPtr->togl_HasAccumulationBuffer) /* 43 */
#endif
#ifndef Togl_HasDestinationAlpha
#define Togl_HasDestinationAlpha \
	(toglStubsPtr->togl_HasDestinationAlpha) /* 44 */
#endif
#ifndef Togl_HasStencilBuffer
#define Togl_HasStencilBuffer \
	(toglStubsPtr->togl_HasStencilBuffer) /* 45 */
#endif
#ifndef Togl_StereoMode
#define Togl_StereoMode \
	(toglStubsPtr->togl_StereoMode) /* 46 */
#endif
#ifndef Togl_HasMultisample
#define Togl_HasMultisample \
	(toglStubsPtr->togl_HasMultisample) /* 47 */
#endif
#ifndef Togl_CopyContext
#define Togl_CopyContext \
	(toglStubsPtr->togl_CopyContext) /* 48 */
#endif

#endif /* defined(USE_TOGL_STUBS) && !defined(USE_TOGL_STUB_PROCS) */

/* !END!: Do not edit above this line. */

#endif
