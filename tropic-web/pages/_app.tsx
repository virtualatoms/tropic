import "@mantine/core/styles.css";
import "@mantine/code-highlight/styles.css";
import "highlight.js/styles/github.css";
import "@/assets/style.css";

import Head from "next/head";
import { MantineProvider } from "@mantine/core";
import {
  CodeHighlightAdapterProvider,
  createHighlightJsAdapter,
} from "@mantine/code-highlight";
import { theme } from "@/theme";
import { Layout } from "@/components/Layout";

import hljs from "highlight.js/lib/core";
import bash from "highlight.js/lib/languages/bash";
import python from "highlight.js/lib/languages/python";

// 2. Register the languages with highlight.js
hljs.registerLanguage("bash", bash);
hljs.registerLanguage("python", python);

const highlightJsAdapter = createHighlightJsAdapter(hljs);

export default function App({ Component, pageProps }: any) {
  return (
    <MantineProvider theme={theme}>
      <CodeHighlightAdapterProvider adapter={highlightJsAdapter}>
        <Head>
          <title>TROPIC</title>
          <meta
            name="viewport"
            content="minimum-scale=1, initial-scale=1, width=device-width, user-scalable=no"
          />
          <link rel="shortcut icon" href="/favicon.svg" />
        </Head>
        <Layout>
          <Component {...pageProps} />
        </Layout>
      </CodeHighlightAdapterProvider>
    </MantineProvider>
  );
}
