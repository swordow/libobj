/*
*	swordow, 2021/01/29
*/
#include <memory>
#include <stdlib.h>
#include "objfile.h"

#define DEBUG_OUTPUT 0
#define SAFE_RELEASE(ptr) if ((ptr)) {free((void*)(ptr));}

#if DEBUG_OUTPUT
#define DebugOutput(fp, fmt, ...) fprintf(fp, fmt, __VA_ARGS__)
#else
#define DebugOutput(fp, fmt, ...)
#endif

typedef struct RangeIter
{
    const char* Begin, *End;
}RangeIter;

typedef struct ObjRawFace
{
    int32_t V[3];
    int32_t T[3];
    int32_t N[3];
}ObjRawFace;

typedef struct ObjRawGroup
{
    int32_t         FCap;
    int32_t         FCount;
    ObjRawFace*     Faces;
    ObjRawGroup*    Next;
    int32_t         SubIdx;
    char            Name[256];
    char            MtlName[256];
}ObjRawGroup;

typedef struct ObjParserState
{
	int32_t		    VCount;
    int32_t         VCap;
	int32_t		    NCount;
    int32_t         NCap;
	int32_t		    TCount;
    int32_t         TCap;
	ObjVector3*     Vertexs;
	ObjVector4*     TexCoords;
	ObjVector3*     Normals;
    const char*     Buffer;
    uint32_t        Bytes;
    uint32_t        Pos;
    char            Error[1024];
    ObjRawGroup*    RawGroups;
    ObjRawGroup*    EndRawGroup;
    FILE*           OutFile;
    uint32_t        LineNum;
    ObjMtlLib       MtlLib;
    ObjMtl*         EndMtl;
}ObjParserState;

struct ObjCommand;

typedef int32_t(*ParserProc)(ObjParserState* state, ObjectFile* f, ObjCommand* cmd, RangeIter r);

typedef struct ObjCommand
{
    char* Cmd;
    uint32_t N;
    ParserProc Proc;
}ObjCommand;

#define OBJCMD(name, proc) {#name, (uint32_t)strlen(#name), proc}

static void* MemAlloc(size_t bytes)
{
    void* b = malloc(bytes);
    if (b == 0)
    {
        printf("Out of memory!");
        exit(1);
    }
    return b;
}

static void MemFree(void* b)
{
    if (b != 0)
    {
        free(b);
    }
}

static void* MemRealloc(void* old, size_t bytes)
{
    void* b = realloc(old, bytes);
    if (b == 0)
    {
        printf("realloc Out of memory!");
        exit(1);
    }
    return b;
}

static void* LoadFileBuffer(const char* file, uint32_t* bytes)
{
    FILE* f = fopen(file, "rb");
    if (!f)
    {
        printf("file %s not found!\n", file);
        return 0;
    }
    uint8_t* buffer = (uint8_t*)MemAlloc(102400u);
    int32_t pos = 0;
    int32_t cacheSize = 1;
    while (1)
    {
        size_t s = fread(&buffer[pos], 1u, 102400u, f);
        if (s >= 102400u)
        {
            cacheSize++;
            pos += (int32_t)s;
            buffer = (uint8_t*)MemRealloc(buffer, cacheSize * 102400u);
            if (!buffer)
            {
                printf("Load file failed! Out of Memory!");
                return 0;
            }
        }
        else
        {
            pos += (int32_t)s;
            break;
        }
    }
    *bytes = pos;
    return buffer;
}

static void FreeFileBuffer(void* buffer)
{
    SAFE_RELEASE(buffer);
}

ObjRawGroup* NewObjRawGroup(ObjParserState* state)
{
    ObjRawGroup* g = (ObjRawGroup*)MemAlloc(sizeof(ObjRawGroup));
    g->Next = 0;
    g->Name[0] = 0;
    g->MtlName[0] = 0;
    g->FCount = 0;
    g->FCap = 1024;
    g->Faces = (ObjRawFace*)MemAlloc(sizeof(ObjRawFace)*g->FCap);
    g->SubIdx = 0;
    return g;
}


static void InitializeParserState(ObjParserState* state, const char* outfileName)
{
    state->VCount       = state->NCount = state->TCount = 0;
    state->VCap         = state->NCap = state->TCap = 1024;
    state->Vertexs      = (ObjVector3*)MemAlloc(sizeof(ObjVector3) * state->VCap);
    state->Normals      = (ObjVector3*)MemAlloc(sizeof(ObjVector3) * state->NCap);
    state->TexCoords    = (ObjVector4*)MemAlloc(sizeof(ObjVector4) * state->TCap);
    state->RawGroups    = 0;
    state->EndRawGroup  = 0;
    state->Error[0]     = 0;
    state->MtlLib.Name[0] = 0;
    state->MtlLib.Mtls = 0;
    state->EndMtl = 0;
#if DEBUG_OUTPUT
    state->OutFile =  fopen(outfileName, "wt");
#else
    state->OutFile = 0;
#endif
}

void FreeObjRawGroup(ObjRawGroup* g)
{
    SAFE_RELEASE(g->Faces);
    SAFE_RELEASE(g);
}

static void FinializeParserState(ObjParserState* state)
{
    state->Buffer = 0;
    state->Bytes = 0;
    state->Pos = 0;
    SAFE_RELEASE(state->Vertexs);
    SAFE_RELEASE(state->TexCoords);
    SAFE_RELEASE(state->Normals);
    ObjRawGroup* g = state->RawGroups;
    while (g)
    {
        ObjRawGroup* next = g->Next;
        FreeObjRawGroup(g);
        g = next;
    }
    state->RawGroups = nullptr;
    state->EndRawGroup = nullptr;
    ObjMtl* mtl = state->MtlLib.Mtls;
    state->MtlLib.Mtls = 0;
    while (mtl)
    {
        ObjMtl* next = mtl->Next;
        SAFE_RELEASE(mtl);
        mtl = next;
    }
    if (state->OutFile)
    {
        fclose(state->OutFile);
        state->OutFile = 0;
    }
}


static void SkipWhitespace(const char* buffer, RangeIter* r)
{
    while (r->Begin < r->End && (*r->Begin == ' ' || *r->Begin == '\t')) r->Begin += 1;
}

static RangeIter ReadLine(const char* buffer, uint32_t* index, uint32_t bytes)
{
    RangeIter r = { &buffer[*index], &buffer[*index] };
    while (*index < bytes && buffer[*index] != '\n' && buffer[*index] != '\r') (*index)++;
    r.End = &buffer[*index];

    if (buffer[*index] == '\r') (*index)++;

    (*index)++;
    return r;
}

static RangeIter ReadData(const char* buffer, RangeIter* r, char sep, int32_t alsoWS)
{
    RangeIter ret = { r->Begin, r->Begin };
    while (r->Begin < r->End)
    {
        char check = *r->Begin;
        if (check == sep || (alsoWS && (check == ' ' || check == '\t')))
        {
            break;
        }
        r->Begin += 1;
    }
    ret.End = r->Begin;
    return ret;
}

ObjMesh* NewObjMesh(ObjParserState* state)
{
    ObjMesh* mesh = (ObjMesh*)MemAlloc(sizeof(ObjMesh));
    mesh->Next = 0;
    mesh->Name[0] = 0;
    mesh->FCount = mesh->VCount = 0;
    return mesh;
}

ObjGroup* NewObjGroup(ObjParserState* state)
{
    ObjGroup* g = (ObjGroup*)MemAlloc(sizeof(ObjGroup));
    g->Next = 0;
    g->Name[0] = 0;
    g->MeshList = NewObjMesh(state);
    return g;
}

ObjMtl* NewObjMtl(ObjParserState* state)
{
    ObjMtl* m = (ObjMtl*)MemAlloc(sizeof(ObjMtl));
    m->MapBump[0] = 0;
    m->MapKd[0] = 0;
    m->Name[0] = 0;
    m->Next = nullptr;
    return m;
}

int32_t Obj_atoi(const char* buffer, RangeIter r)
{
    char num[128] = { 0 };
    strncpy(num, r.Begin, size_t(r.End - r.Begin));
    return atoi(num);
}

float Obj_atof(const char* buffer, RangeIter r)
{
    char num[128] = { 0 };
    strncpy(num, r.Begin, size_t(r.End - r.Begin));
    return (float)atof(num);
}

void ObjRawGroupAddFace(ObjRawGroup* group, ObjRawFace face)
{
    if (group->FCount + 1 >= group->FCap)
    {
        group->Faces = (ObjRawFace*)MemRealloc(group->Faces, (group->FCap + 1024) * sizeof(ObjRawFace));
        group->FCap += 1024;
    }
    group->Faces[group->FCount] = face;
    group->FCount += 1;
}

static void ParseState_AddVertex(ObjParserState* state, ObjVector3* v)
{
    if (state->VCount + 1 == state->VCap)
    {
        state->Vertexs = (ObjVector3*)MemRealloc(state->Vertexs, (state->VCap + 1024)*sizeof(ObjVector3));
        state->VCap += 1024;
    }
    state->Vertexs[state->VCount] = *v;
    state->VCount += 1;
}

static void ParseState_AddNormal(ObjParserState* state, ObjVector3* v)
{
    if (state->NCount + 1 == state->NCap)
    {
        state->Normals = (ObjVector3*)MemRealloc(state->Normals, (state->NCap + 1024) * sizeof(ObjVector3));
        state->NCap += 1024;
    }
    state->Normals[state->NCount] = *v;
    state->NCount += 1;
}

static void ParseState_AddTexcoord(ObjParserState* state, ObjVector4* v)
{
    if (state->TCount + 1 == state->TCap)
    {
        state->TexCoords = (ObjVector4*)MemRealloc(state->TexCoords, (state->TCap + 1024) * sizeof(ObjVector4));
        state->TCap += 1024;
    }
    state->TexCoords[state->TCount] = *v;
    state->TCount += 1;
}

static void ParseState_ObjRawGroup(ObjParserState* state, ObjRawGroup* g)
{
    if (state->RawGroups == 0)
    {
        state->RawGroups = g;
        state->EndRawGroup = g;
    }
    else
    {
        state->EndRawGroup->Next = g;
        state->EndRawGroup = g;
    }
}

int32_t Parse_V(ObjParserState* state, ObjectFile* f, ObjCommand* cmd, RangeIter r)
{
    int i = 0;
    ObjVector3 point = { 0,0,0 };
    // skip 'v'
    r.Begin += cmd->N;
    for (i = 0; i < 3; i++)
    {
        SkipWhitespace(state->Buffer, &r);
        RangeIter range = ReadData(state->Buffer, &r, ' ', true);
        if (range.Begin >= range.End)
        {
            sprintf(state->Error, "Line:%d - v has no data!", state->LineNum);
        }
        point.v[i] = Obj_atof(state->Buffer, range);
    }
    ParseState_AddVertex(state, &point);
    DebugOutput(state->OutFile, "v  %0.4f %0.4f %0.4f\n", point.x, point.y, point.z);
    return 0;
}

int32_t Parse_VN(ObjParserState* state, ObjectFile* f, ObjCommand* cmd, RangeIter r)
{
    int i = 0;
    ObjVector3 point = { 0,0,0 };
    // skip 'vn'
    r.Begin += cmd->N;
    for (i = 0; i < 3; i++)
    {
        SkipWhitespace(state->Buffer, &r);
        RangeIter range = ReadData(state->Buffer, &r, ' ', true);
        if (range.Begin >= range.End)
        {
            sprintf(state->Error, "Line:%d - vn has no data!", state->LineNum);
        }
        point.v[i] = Obj_atof(state->Buffer, range);
    }
    ParseState_AddNormal(state, &point);
    DebugOutput(state->OutFile, "vn %0.4f %0.4f %0.4f\n", point.x, point.y, point.z);
    return 0;
}

int32_t Parse_VT(ObjParserState* state, ObjectFile* f, ObjCommand* cmd, RangeIter r)
{
    int i = 0;
    ObjVector4 point = { 0,0,0 };
    // skip 'vn'
    r.Begin += cmd->N;
    for (i = 0; i < 2; i++)
    {
        char num[128] = { 0 };
        // skip
        SkipWhitespace(state->Buffer, &r);
        RangeIter range = ReadData(state->Buffer, &r, ' ', true);
        if (range.Begin >= range.End)
        {
            sprintf(state->Error, "Line:%d - vt has no data!", state->LineNum);
        }
        point.v[i] = Obj_atof(state->Buffer, range);
    }

    // try if it has w
    SkipWhitespace(state->Buffer, &r);
    RangeIter range = ReadData(state->Buffer, &r, ' ', true);
    if (range.Begin < range.End)
    {
        point.v[2] = Obj_atof(state->Buffer, range);
    }

    ParseState_AddTexcoord(state, &point);
    DebugOutput(state->OutFile, "vt %0.4f %0.4f %0.4f\n", point.x, point.y, point.z);
    return 0;
}
int32_t Parse_SubF(ObjParserState* state, ObjectFile* f, RangeIter* r, ObjRawFace* face, int32_t fItem)
{
    int32_t mode = 0x01;

    /*
        v
        v/vt
        v/vt/vn
        v//vn
    */
    SkipWhitespace(state->Buffer, r);
    RangeIter range = ReadData(state->Buffer, r, '/', true);
    // must have vertex index
    if (range.Begin >= range.End)
    {
        return 0;
    }

    face->V[fItem] = Obj_atoi(state->Buffer, range);

    // check has vt
    if (*r->Begin == '/')
    {
        // skip '/'
        r->Begin += 1;
        range = ReadData(state->Buffer, r, '/', true);
        if (range.Begin < range.End)
        {
            face->T[fItem] = Obj_atoi(state->Buffer, range);
            mode |= 0x02;
        }

        // check has vn 
        if (*r->Begin == '/')
        {
            // skip '/'
            r->Begin += 1;
            range = ReadData(state->Buffer, r, ' ', true);
            if (range.Begin < range.End)
            {
                face->N[fItem] = Obj_atoi(state->Buffer, range);
                mode |= 0x04;
            }
            // v//?
            else
            {
                return 0;
            }
        }
        // v/?
        else if ((mode & 0x02) == 0)
        {
            return 0;
        }
    }
    return mode;
}

int32_t Parse_F(ObjParserState* state, ObjectFile* f, ObjCommand* cmd, RangeIter r)
{
    ObjRawFace face = {
        {0,0,0},
        {0,0,0},
        {0,0,0}
    };
   
    // skip 'f'
    r.Begin += cmd->N;

    int32_t mode = Parse_SubF(state, f, &r, &face, 0);
    if (mode == 0)
    {
        sprintf(state->Error, "Line:%d - f format error", state->LineNum);
        return -1;
    }

    if (Parse_SubF(state, f, &r, &face, 1) != mode)
    {
        sprintf(state->Error, "Line:%d - f format not consists!", state->LineNum);
        return -1;
    }
    if (Parse_SubF(state, f, &r, &face, 2) != mode)
    {
        sprintf(state->Error, "Line:%d - f format not consists!", state->LineNum);
        return -1;
    }
    ObjRawGroupAddFace(state->EndRawGroup, face);
    DebugOutput(state->OutFile, "f %d/%d/%d %d/%d/%d %d/%d/%d \n", face.V[0], face.T[0], face.N[0], face.V[1], face.T[1], face.N[1], face.V[2], face.T[2], face.N[2]);
    return 0;
}

int32_t Parse_G(ObjParserState* state, ObjectFile* f, ObjCommand* cmd, RangeIter r)
{
    ObjRawGroup* g = NewObjRawGroup(state);

    // skip 'g'
    r.Begin += cmd->N;
    SkipWhitespace(state->Buffer, &r);
    RangeIter range = ReadData(state->Buffer, &r, ' ', true);
    if (range.Begin >= range.End)
    {
        sprintf(state->Error, "Line:%d - g has no name!", state->LineNum);
        return -1;
    }
    strncpy(g->Name, range.Begin, size_t(range.End - range.Begin));
    g->Name[range.End - range.Begin] = 0;
    ParseState_ObjRawGroup(state, g);
    DebugOutput(state->OutFile, "g %s\n", state->EndRawGroup->Name);
    return 0;
}

int32_t Parse_UseMtl(ObjParserState* state, ObjectFile* f, ObjCommand* cmd, RangeIter r)
{
    // skip 'usemtl'
    r.Begin += cmd->N;

    SkipWhitespace(state->Buffer, &r);
    RangeIter range = ReadData(state->Buffer, &r, ' ', true);
    
    if (range.Begin >= range.End)
    {
        sprintf(state->Error, "Line:%d - usemtl has no name!", state->LineNum);
        return -1;
    }

    // crate new group
    if (state->EndRawGroup->MtlName[0] != 0)
    {
        ObjRawGroup* ng = NewObjRawGroup(state);
        strcpy(ng->Name, state->EndRawGroup->Name);
        sprintf(&ng->Name[strlen(ng->Name)], "_%d", state->EndRawGroup->SubIdx + 1);
        ng->SubIdx = state->EndRawGroup->SubIdx + 1;
        ParseState_ObjRawGroup(state, ng);
    }

    strncpy(state->EndRawGroup->MtlName, range.Begin, size_t(range.End - range.Begin));
    state->EndRawGroup->MtlName[range.End - range.Begin] = 0;
    DebugOutput(state->OutFile, "usemtl %s\n", state->EndRawGroup->MtlName);
    return 0;
}

int32_t Parse_NewMtl(ObjParserState* state, ObjectFile* f, ObjCommand* cmd, RangeIter r)
{
    // skip 'newmtl'
    r.Begin += cmd->N;

    SkipWhitespace(state->Buffer, &r);
    RangeIter range = ReadData(state->Buffer, &r, ' ', true);

    if (range.Begin >= range.End)
    {
        sprintf(state->Error, "Line:%d - newmtl has no name!", state->LineNum);
        return -1;
    }

    ObjMtl* mtl = NewObjMtl(state);
    strncpy(mtl->Name, range.Begin, size_t(range.End - range.Begin));
    mtl->Name[range.End - range.Begin] = 0;

    if (state->EndMtl == 0)
    {
        state->EndMtl = mtl;
        state->MtlLib.Mtls = mtl;
    }
    else
    {
        state->EndMtl->Next = mtl;
        state->EndMtl = mtl;
    }
    DebugOutput(state->OutFile, "newmtl %s\n", mtl->Name);
    return 0;
}

int32_t Parse_MtlMap(ObjParserState* state, ObjectFile* f, ObjCommand* cmd, RangeIter r)
{
    // ignore all options...
    // only need filename
    // skip 'newmtl'
    r.Begin += cmd->N;

    SkipWhitespace(state->Buffer, &r);
    RangeIter range = ReadData(state->Buffer, &r, ' ', true);

    if (range.Begin >= range.End)
    {
        sprintf(state->Error, "Line:%d - %s has no name!",state->LineNum, cmd->Cmd);
        return -1;
    }

    if (strcmp(cmd->Cmd, "map_Kd") == 0)
    {
        strncpy(state->EndMtl->MapKd, range.Begin, size_t(range.End - range.Begin));
        state->EndMtl->MapKd[range.End - range.Begin] = 0;
        DebugOutput(state->OutFile, "%s %s\n", cmd->Cmd, state->EndMtl->MapKd);
    }
    else if(strcmp(cmd->Cmd, "map_bump") == 0)
    {
        strncpy(state->EndMtl->MapBump, range.Begin, size_t(range.End - range.Begin));
        state->EndMtl->MapBump[range.End - range.Begin] = 0;
        DebugOutput(state->OutFile, "%s %s\n", cmd->Cmd, state->EndMtl->MapBump);
    }
    return 0;
}

struct NamedVertex
{
    char FName[256];
    ObjVector3 V;
    ObjVector4 T;
    ObjVector3 N;
};


int32_t GenerateObjectFile(ObjParserState* state, ObjectFile* of)
{
    ObjGroup* tmpRoot = NewObjGroup(state);
    ObjGroup* curGroup = tmpRoot;
    ObjMesh* curMesh = 0;
    ObjRawGroup* rg = state->RawGroups;
    while (rg)
    {
        if (rg->SubIdx == 0)
        {
            curGroup->Next = NewObjGroup(state);
            curGroup = curGroup->Next;
            curMesh = curGroup->MeshList;
            strcpy(curGroup->Name, rg->Name);
            strcpy(curMesh->Name, rg->Name);
        }
        else
        {
            curMesh->Next = NewObjMesh(state);
            curMesh = curMesh->Next;
            strcpy(curMesh->Name, rg->Name);
        }

        ObjMtl* mtl = state->MtlLib.Mtls;
        while (mtl)
        {
            if (strcmp(mtl->Name, rg->MtlName) == 0)
            {
                break;
            }
            mtl = mtl->Next;
        }
        curMesh->Mtl = mtl;
        
        curMesh->FCount = rg->FCount;
        curMesh->Faces = (ObjFace*)MemAlloc(sizeof(ObjFace) * rg->FCount);
        int32_t nvCount = 0;
        int32_t nvCap = 1024;
        NamedVertex* nv = (NamedVertex*)MemAlloc(sizeof(NamedVertex) * 1024);
        for (int32_t i = 0; i < rg->FCount; i++)
        {
            for (int32_t j = 0; j < 3; j++)
            {
                char Name[256] = { 0 };
                sprintf(Name, "%d_%d_%d", rg->Faces[i].V[j], rg->Faces[i].T[j], rg->Faces[i].N[j]);

                // find 
                int32_t n = 0;
                for (n = 0; n < nvCount; n++)
                {
                    // found it, break
                    if (strcmp(Name, nv[n].FName) == 0)
                    {
                        break;
                    }
                }

                // found
                if (n < nvCount)
                {
                    curMesh->Faces[i].P[j] = n;
                }
                // not found, invert new point to nv
                else
                {
                    // inser to nv
                    curMesh->Faces[i].P[j] = nvCount;

                    if (nvCount + 1 >= nvCap)
                    {
                        nv = (NamedVertex*)MemRealloc(nv, sizeof(NamedVertex)*(nvCap + 1024));
                        nvCap += 1024;
                    }

                    strcpy(nv[nvCount].FName, Name);
                    nv[nvCount].V = state->Vertexs[rg->Faces[i].V[j]-1];
                    nv[nvCount].T = state->TexCoords[rg->Faces[i].T[j]-1];
                    nv[nvCount].N = state->Normals[rg->Faces[i].N[j]-1];
                    nvCount += 1;
                }
            }
            
            // compute face normal
            ObjVector3 N0 = state->Normals[rg->Faces[i].N[0] - 1];
            ObjVector3 N1 = state->Normals[rg->Faces[i].N[1] - 1];
            ObjVector3 N2 = state->Normals[rg->Faces[i].N[2] - 1];

            curMesh->Faces[i].N = ObjVector3{
                (N0.x + N1.x + N2.x) / 3.0f,
                (N0.y + N1.y + N2.y) / 3.0f,
                (N0.z + N1.z + N2.z) / 3.0f
            };
        }
        curMesh->VCount = nvCount;
        curMesh->UVLayerCount = 1;
        curMesh->Points = (ObjVector3*)MemAlloc(sizeof(ObjVector3) * nvCount);
        curMesh->UVs[0] = (ObjVector4*)MemAlloc(sizeof(ObjVector4) * nvCount);
        for (int32_t i = 0; i < nvCount; i++)
        {
            curMesh->Points[i] = nv[i].V;
            curMesh->UVs[0][i] = nv[i].T;
        }
        MemFree(nv);
        rg = rg->Next;
    }
    of->GroupList = tmpRoot->Next;
    of->MtlLib = state->MtlLib;
    state->MtlLib.Mtls = 0;
    state->MtlLib.Name[0] = 0;
    return 0;
}

int32_t ProcessCommands(ObjParserState* state, ObjectFile* of, const ObjCommand* objcmds)
{
    state->LineNum = 1;
    RangeIter line = ReadLine(state->Buffer, &state->Pos, state->Bytes);
    for (; line.End < (state->Buffer + state->Bytes); line = ReadLine(state->Buffer, &state->Pos, state->Bytes))
    {
        int32_t cmdIdx = 0;
        ObjCommand cmd = objcmds[cmdIdx];
        int32_t processed = 0;
        
        if (*line.Begin == ' ' || *line.Begin == '\t')
        {
            SkipWhitespace(state->Buffer, &line);
        }

        while (cmd.Cmd != 0)
        {
            if (strncmp(cmd.Cmd, line.Begin, cmd.N) == 0)
            {
                if (cmd.Proc)
                {
                    if (cmd.Proc(state, of, &cmd, line) != 0)
                    {
                        printf("Error %s\n", state->Error);
                        return -1;
                    }
                    else
                    {
                        processed = 1;
                        break;
                    }
                }
            }
            cmdIdx++;
            cmd = objcmds[cmdIdx];
        }

        if (!processed)
        {
            char ll[1024] = { 0 };
            strncpy(ll, line.Begin, size_t(line.End - line.Begin));
            DebugOutput(state->OutFile, "%s\n", ll);
        }
        state->LineNum += 1;
    }
    return 0;
}

const static ObjCommand MtlCmds[] =
{
    OBJCMD(newmtl, &Parse_NewMtl),
    OBJCMD(map_Kd, &Parse_MtlMap),
    OBJCMD(map_bump, &Parse_MtlMap),
    {0, 0, 0}
};

int32_t Parse_MtlLib(ObjParserState* state, ObjectFile* f, ObjCommand* cmd, RangeIter r)
{
    // skip 'usemtl'
    r.Begin += cmd->N;
    SkipWhitespace(state->Buffer, &r);
    RangeIter range = ReadData(state->Buffer, &r, ' ', true);
    if (range.Begin >= range.End)
    {
        sprintf(state->Error, "Line:%d - mtlib has no name!", state->LineNum);
        return -1;
    }
    if (state->MtlLib.Name[0] != 0)
    {
        sprintf(state->Error, "Line:%d - only supports one mtlib!", state->LineNum);
        return -1;
    }
    strncpy(state->MtlLib.Name, range.Begin, size_t(range.End - range.Begin));
    state->MtlLib.Name[range.End - range.Begin] = 0;
    DebugOutput(state->OutFile, "mtllib %s\n", state->MtlLib.Name);

    char mtlFile[256] = { 0 };
    const char* lastSlash = strrchr(f->File, '/');
    strncpy(mtlFile, f->File, size_t(lastSlash-f->File)+1);
    strcat(mtlFile, state->MtlLib.Name);
    printf("Load MTLLib %s\n", mtlFile);
    uint32_t mtlBytes = 0;
    void* mtlBuffer = LoadFileBuffer(mtlFile, &mtlBytes);
    if (!mtlBuffer)
    {
        sprintf(state->Error, "Line:%d - mtllib %s not found!", state->LineNum, mtlFile);
        return -1;
    }
    ObjParserState mtlState;
    InitializeParserState(&mtlState, "objmtl.mtl");
    mtlState.Buffer = (const char*)mtlBuffer;
    mtlState.Bytes = mtlBytes;
    mtlState.Pos = 0;
    int32_t ret = ProcessCommands(&mtlState, nullptr, MtlCmds);
    if (ret != 0)
    {
        sprintf(state->Error, "Parse Mtllib failed!\n");
        FinializeParserState(&mtlState);
        FreeFileBuffer(mtlBuffer);
        return -1;
    }

    state->EndMtl = mtlState.EndMtl;
    state->MtlLib = mtlState.MtlLib;

    mtlState.EndMtl = 0;
    mtlState.MtlLib.Mtls = 0;
    FinializeParserState(&mtlState);
    FreeFileBuffer(mtlBuffer);
    return 0;
}

const static ObjCommand ObjCmds[] =
{
    OBJCMD(usemtl, &Parse_UseMtl),
    OBJCMD(mtllib, &Parse_MtlLib),
    OBJCMD(vt, &Parse_VT),
    OBJCMD(vn, &Parse_VN),
    OBJCMD(v, &Parse_V),
    OBJCMD(f, &Parse_F),
    OBJCMD(s, 0),
    OBJCMD(g, &Parse_G),
    {0, 0, 0}
};

int32_t ParseObjectFile(ObjectFile* of, const char* buffer, int32_t bytes)
{
    ObjParserState state;

    InitializeParserState(&state, "objfile.obj");
    state.Buffer = buffer;
    state.Pos = 0;
    state.Bytes = bytes;
   
    int32_t ret = ProcessCommands(&state, of, ObjCmds);

    if (ret != 0)
    {
        return -1;
    }

    GenerateObjectFile(&state, of);
    FinializeParserState(&state);
    return 0;
}

ObjectFile* LoadObjectFile(const char* file)
{
    ObjectFile* of = (ObjectFile*)MemAlloc(sizeof(ObjectFile));
    uint32_t bytes = 0;
    void* buffer = LoadFileBuffer(file, &bytes);
    if (buffer == 0)
    {
        printf("Load file failed!\n");
        MemFree(of);
        return 0;
    }
    strcpy(of->File, file);
    if (ParseObjectFile(of, (const char*)buffer, bytes) != 0)
    {
        printf("Parse file failed!\n");
        FreeFileBuffer(buffer);
        MemFree(of);
        return 0;
    }
    FreeFileBuffer(buffer);
    return of;
}
