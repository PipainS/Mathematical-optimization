using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml;

namespace OROM
{
    public class OptimalResult
    {
        public Vector xopt { get; set; }
        public Vector ogr { get; set; }
        public Vector fopt { get; set; }
        public int CountIterations = 0;
        public OptimalResult(Vector x, Vector rest, Vector f, int CountIter = 0)
        {
            xopt = x.Copy();
            ogr = null;
            if (rest != null) ogr = rest.Copy();
            fopt = f.Copy();
            CountIterations = CountIter;
        }
        public OptimalResult()
        {
            xopt = null;
            ogr = null;
            fopt = null;
            CountIterations = 0;
        }
        public void ToXML(string NameFile)
        {
            XmlWriterSettings settings = new XmlWriterSettings();
            settings.Indent = true;
            settings.IndentChars = (" ");
            XmlWriter writer = XmlWriter.Create(NameFile, settings);
            writer.WriteStartElement("OptimumResult");
            for (int i = 0; i < xopt.Size; i++)
            {
                writer.WriteStartElement("Xopt");
                writer.WriteAttributeString("Index", i.ToString());
                writer.WriteString(xopt[i].ToString());
                writer.WriteEndElement();
            }
            if (ogr != null)
            {
                for (int i = 0; i < ogr.Size; i++)
                {
                    writer.WriteStartElement("Ogr");
                    writer.WriteAttributeString("Index", i.ToString());
                    writer.WriteString(ogr[i].ToString());
                    writer.WriteEndElement();
                }
            }
            if (fopt.Size == 1)
            {
                writer.WriteElementString("OptimumFunction", fopt[0].ToString());
            }
            else
            {
                for (int i = 0; i < fopt.Size; i++)
                {
                    writer.WriteStartElement("F");
                    writer.WriteAttributeString("Index", i.ToString());
                    writer.WriteString(fopt[i].ToString());
                    writer.WriteEndElement();
                }
            }
            writer.WriteEndElement();
            writer.Flush();
            writer.Close();
        }

        public override string ToString()
        {
            if (ogr == null) return string.Format(" Result xopt={0}  fopt={1} count={2}", xopt, fopt, CountIterations);
            return string.Format(" Result xopt={0}  ogr={1}  fopt={2} count={3}", xopt, ogr, fopt, CountIterations);
        }
    }
}
