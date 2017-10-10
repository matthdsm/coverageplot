package coverageplot

import (
	"encoding/json"
	"github.com/gorilla/mux"
	"github.com/meatballhat/negroni-logrus"
	"github.com/urfave/negroni"
	"net/http"
	"strconv"
	"strings"
)

func Serve() *http.Server {

	srv := &http.Server{Addr: ":8080"}
	r := mux.NewRouter()

	r.HandleFunc("/regions", getRegionsHandler)
	r.PathPrefix("/").Handler(http.FileServer(http.Dir("public/")))

	n := negroni.New()
	n.Use(negronilogrus.NewMiddleware())
	n.UseHandler(r)

	http.Handle("/", n)

	//go func() {
	//	if err := srv.ListenAndServe(); err != nil {
	//		// cannot panic, because this probably is an intentional close
	//		log.Printf("Httpserver: ListenAndServe() error: %s", err)
	//	}
	//}()

	return srv
}

func getRegionsHandler(w http.ResponseWriter, r *http.Request) {
	w.Header().Set("Content-Type", "application/json")
	var intervals []BedInterval
	var coverageBed string
	var intervalsFile string
	intervalsFile = r.FormValue("intervalsFile")
	coverageBed = r.FormValue("coverageBed")
	regionParam := r.FormValue("region")

	if regionParam != "" && intervalsFile == "" {
		toks := strings.FieldsFunc(regionParam, SplitRegion)
		s, _ := strconv.Atoi(toks[1])
		e, _ := strconv.Atoi(toks[2])
		info := []BedInfo{{"adhoc", "adhoc", "adhoc", "adhoc", "adhoc"}}
		intervals = append(intervals, BedInterval{toks[0], s, e, info})
		intervalsFile = "intervals"
	}
	if len(intervals) < 1 && intervalsFile != "" {
		intervals = ParseRegions(intervalsFile)
	}
	regions := CalculateStatistics(intervals, intervalsFile, coverageBed, "adhoc")
	j, _ := json.Marshal(regions)
	w.Write(j)
}

func SplitRegion(r rune) bool {
	return r == ':' || r == '-'
}
