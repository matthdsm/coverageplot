package coverageplot

import (
	"encoding/json"
	"github.com/gorilla/mux"
	"github.com/meatballhat/negroni-logrus"
	"github.com/urfave/negroni"
	"log"
	"net/http"
)

func Serve() *http.Server {

	srv := &http.Server{Addr: ":8080"}
	r := mux.NewRouter()

	r.HandleFunc("/regions", getRegionsHandler)
	r.PathPrefix("/").Handler(http.FileServer(http.Dir("../public/")))

	n := negroni.New()
	n.Use(negronilogrus.NewMiddleware())
	n.UseHandler(r)

	http.Handle("/", n)

	go func() {
		if err := srv.ListenAndServe(); err != nil {
			// cannot panic, because this probably is an intentional close
			log.Printf("Httpserver: ListenAndServe() error: %s", err)
		}
	}()

	return srv
}

func getRegionsHandler(w http.ResponseWriter, r *http.Request) {
	w.Header().Set("Content-Type", "application/json")
	var intervalsFile string
	var coverageBed string
	var panelName string
	regions := CalculateStatistics(intervalsFile, coverageBed, panelName)
	j, _ := json.Marshal(regions)
	w.Write(j)
}
