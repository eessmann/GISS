import sys
import os
import os.path
import re
import tempfile

def generated(lines):
    for line in lines:
        record = line.split()
        if len(record) > 3 and \
               record[0] == "#" and record[1] == "Generated" and record[2] == "files:":
            return record[3:]
    return []

def dictionary(d,file):
    p = re.compile(r">\w*_class</A\n")
    instruct = 0
    for line in file:
        record = line.split()
        if re.match(p,line) != None:
            klass = ""
            cap = 1
            for c in line[1:-9]:
                if c == "_":
                    cap = 1
                elif cap:
                    klass += c.capitalize()
                    cap = 0
                else:
                    klass += c
            if klass == "GfsOutputPpm":
                klass = "GfsOutputPPM"
            d[klass] = href
        elif len(record) == 2 and record[0] == "struct" and record[1] == "<A":
            instruct = 1
        elif len(record) == 1 and record[0] == ">;":
            instruct = 0
        elif instruct:
            if line[0:5] == "HREF=":
                val = line[6:-2]
            elif line[0] == ">" and line[-4:] == "</A\n":
                d[line[1:-4]] = val
        elif line[0:5] == "HREF=":
            href = line[6:line.find("#")]

class Example:
    def __init__(self,path):
        if path[0:2] == "./":
            path = path[2:]
        self.path, self.name = os.path.split(path)
        if self.name == "":
            self.name = self.path
        elif self.path == "":
            self.path = self.name
        else:
            self.path += "/" + self.name
        self.section = ["\\subsection","\\subsubsection"][self.path.count("/")]
        file = open(self.path + "/" + self.name + ".gfs")
        lines = file.readlines()
        self.generated = generated(lines)
        if os.access(self.path + "/status", os.R_OK):
            self.status = open(self.path + "/status").readline()
            self.generated.append("status")
        else:
            self.status = None
        p = re.compile(r"\\label\{[a-zA-Z0-9_\-]*\}")
        labels = []
        for line in lines:
            for l in re.findall(p,line):
                labels.append(l[7:-1])

        # adds the full path to references to generated files and makes labels absolute
        lines1 = []
        path = self.path[2:].replace("/", "-")
        for line in lines:
            for gen in self.generated:
                line = line.replace("{" + gen + "}", "{" + self.path + "/" + gen + "}")
            for l in labels:
                line = line.replace("{" + l + "}", "{" + path + "-" + l + "}")
            lines1.append(line)
        lines = lines1

        self.title = []
        self.description = []

        insthg = None
        for line in lines:
            record = line.split()
            if len(record) > 0 and record[0] == "#":
                if len(record) > 1:
                    if record[1] == "Title:":
                        self.title.append(" ".join(record[2:]))
                        insthg = self.title
                    elif record[1] == "Description:":
                        insthg = self.description
                    elif record[1] == "Required" and record[2] == "files:":
                        self.required = record[3:]
                        insthg = None
                    elif record[1] == "Command:":
                        self.command = " ".join(record[2:])
                        insthg = None
                    elif record[1] == "Author:":
                        self.author = " ".join(record[2:])
                        insthg = None
                    elif record[1] == "Running" and record[2] == "time:":
                        self.time = " ".join(record[3:])
                        insthg = None
                    elif record[1] == "Version:":
                        self.version = " ".join(record[2:])
                        insthg = None
                    elif not insthg == None:
                        insthg.append(" ".join(record[1:]))
                elif not insthg == None:
                    insthg.append(" ".join(record[1:]))

        if os.access(self.path + "/runtime", os.R_OK):
            self.runtime = float(open(self.path + "/runtime").readline())
            self.time = ""
            m = int(self.runtime/60.)
            if m > 0:
                self.time += repr(m) + " minutes"
            s = int(self.runtime-60.*m)
            if s > 0:
                self.time += " " + repr(s) + " seconds"
            self.generated.append("runtime")
        else:
            self.runtime = None
            
    def write(self,dico,file=None):
        if file == None:
            file = open(self.path + "/" + self.name + ".tex", 'w')
        if self.status:
            file.write(self.section + "{" + self.status + "\n".join(self.title) + "}\n")
        else:
            file.write(self.section + "{" + "\n".join(self.title) + "}\n")
        file.write("\\begin{description}\n")
        file.write("\\item[Author]" + self.author + "\n")
        file.write("\\item[Command]" + "{\\tt " + self.command + "}\n")
        file.write("\\item[Version]" + self.version + "\n")
        f = self.name + ".gfs"
        required = " " + f + \
                   " \\htmladdnormallinkfoot{(view)}{" + self.path + "/" + f + ".html}" +\
                   " \\htmladdnormallinkfoot{(download)}{" + self.path + "/" + f + "}\\\\"
        for f in self.required:
            required += " \\htmladdnormallinkfoot{" + f + "}{" + self.path + "/" + f + "}"
        file.write("\\item[Required files]" + required + "\n")
        file.write("\\item[Running time]" + self.time + "\n")
        file.write("\\end{description}\n")
        file.write("\n".join(self.description))
        self.colorize(dico)

    def colorize(self,dico):
        file = open(self.path + "/" + self.name + ".gfs")
        out = open(self.path + "/" + self.name + ".gfs.html", 'w')
        path = "../" * (self.path.count("/") + 3) + "reference/"
        out.write("<html><head><title>\n" + self.name + ".gfs")
        out.write("</title></head><body><tt>\n")
        infile = insthg = 0
        for line in file:
            l = ""
            first = 1
            comment = 0
            for c in line:
                if c == " " and first:
                    l += "&nbsp; "
                else:
                    if first and c == "#":
                        comment = 1
                    first = 0
                    if c == '<':
                        l += "&lt"
                    elif c == '>':
                        l += "&gt"
                    elif c == '&':
                        l += "&amp"
                    else:
                        l += c
            if comment:
                if infile:
                    record = line.split()
                    if insthg or len(record) > 1:
                        out.write("<font color=\"#5285B4\">")
                        out.write(l)
                        out.write("</font><br>")
                        insthg = 1
                else:
                    record = line.split()
                    if len(record) > 2 and record[1] == "Generated" and record[2] == "files:":
                        infile = 1
                        insthg = 0
            else:
                record = l.split()
                for r in record:
                    key = None
                    if dico.has_key(r):
                        key = r
                    elif dico.has_key("Gfs" + r):
                        key = "Gfs" + r
                    if key != None:
                        out.write("<a href=\"" + path + dico[key] + "\">" + r + "</a> ")
                    else:
                        out.write(r + " ")
                out.write("<br>\n")
        out.write("</tt></body>\n</html>\n")

    def test(self):
        wdname = tempfile.mkdtemp()
        path = os.getcwd() + "/" + self.path + "/"
        files = path + self.name + ".gfs"
        for f in self.required:
            files += " " + path + f
        out = os.popen("cd " + wdname + " && " +\
                       "mkdir test && cd test && " +\
                       "cp -f " + files + " . && " +\
                       "awk '{ if ($1 == \"Time\" || $1 == \"GfsTime\")" +\
                       "  print $0 \"\\nTime { iend = 1 }\";" +
                       "else print $0;"
                       "}' < " + self.name + ".gfs > " + self.name + ".tmp && " +\
                       "mv -f " + self.name + ".tmp " + self.name + ".gfs && ( " +\
                       self.command + " ) 2>&1")
        lines = out.readlines()
        status = out.close()
        os.system("rm -r -f " + wdname)
        if status != None:
            return status,lines
        else:
            return None,None

    def run(self,env=""):
        out = os.popen("cd " + self.path + " && " +\
                       "sh -c \"time -p " + env + " " + self.command + "\" 2>&1")
        lines = []
        for l in out:
            record = l.split()
            if len(record) > 0:
                if record[0] == "user":
                    self.runtime = float(record[1])
                    print >>open(self.path + "/runtime",'w'), self.runtime
                elif record[0] != "real" and record[0] != "sys":
                    lines.append(l)
            else:
                lines.append(l)
        status = out.close()
        if status != None:
            return status,lines
        else:
            return None,None
