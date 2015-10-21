#ifndef MTRACE_H
#define MTRACE_H

#include <cstring>
#include <stdio.h>
#include <utility>

#include <ai.h>

typedef void (*MsgFunc) (const char* format, ...);

template<typename ...Args>
void MsgFuncTrace(const MsgFunc msgFunc,
                  const char* file,
                  const int line,
                  const char* func,
                  const AtNode* node,
                  const char* format,
                  Args &&...args)
{
    char prefixedFormat[BUFSIZ] = {0};
    strcat (prefixedFormat, "[%s] "); // type
    strcat (prefixedFormat, "[%s | %d | %s] "); // file, line, and func
    strcat (prefixedFormat, format);
    (*msgFunc)(prefixedFormat, AiNodeEntryGetName(AiNodeGetNodeEntry(node)), file, line, func, std::forward<Args>(args)...);
}

template<typename ...Args>
void MsgFuncPrefix(const MsgFunc msgFunc,
                   const AtNode* node,
                  const char* format,
                   Args &&...args)
{
    char prefixedFormat[BUFSIZ] = {0};
    strcat (prefixedFormat, "[%s] "); // type
    strcat (prefixedFormat, format);
    (*msgFunc)(prefixedFormat, AiNodeEntryGetName(AiNodeGetNodeEntry(node)), std::forward<Args>(args)...);
}

#ifdef NDEBUG
#define MTRACE(FUNC, NODE, FORMAT, ...) MsgFuncPrefix (FUNC, NODE, FORMAT, ##__VA_ARGS__)
#else
#define MTRACE(FUNC, NODE, FORMAT, ...) MsgFuncTrace (FUNC, __FILE__, __LINE__, __PRETTY_FUNCTION__, NODE, FORMAT, ##__VA_ARGS__)
#endif

#endif // MTRACE_H
